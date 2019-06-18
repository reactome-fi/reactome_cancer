/*
 * Created on Nov 19, 2011
 *
 */
package org.reactome.cancer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.fi.MCLClusteringHelper;
import org.reactome.r3.fi.MCLClusteringResult;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * @author gwu
 *
 */
public class MCLClusterWrapper {
    private MCLClusteringHelper mclHelper;
    private double inflation = 5.0d;
    private boolean useAbsolute = true;
    
    public MCLClusterWrapper() {
        mclHelper = new MCLClusteringHelper();
        mclHelper.setMclScript(R3Constants.mclScript);
        mclHelper.setTempDirName(R3Constants.TEMP_DIR);
    }
    
    public void setTempDirName(String dirName) {
        mclHelper.setTempDirName(dirName);
    }
    
    public void setInflation(double inflation) {
        this.inflation = inflation;
    }
    
    public void setUseAbsolute(boolean useAbsolute) {
        this.useAbsolute = useAbsolute;
    }
    
    public void setKeepTempFile(boolean keep) {
        mclHelper.setKeepTempFile(keep);
    }
    
    public List<Set<String>> mclCluster(String geneExpFileName,
                                        int sizeCutoff,
                                        double corCutoff) throws Exception {
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(geneExpFileName);
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fisWithCorrs = helper.calculateGeneExpCorrForFIs(geneToSampleToValue,
                                                                     fis, 
                                                                     this.useAbsolute,
                                                                     null);
        return mclCluster(fisWithCorrs, 
                          fis,
                          sizeCutoff,
                          corCutoff);
    }
    
    /**
     * Do a MCL clustering by using p-values for genes. Gene-based p-values are converted
     * into edge-weights by using geometric or arithemtic mean of two genes involved in the same FIs. 
     * If p-values for either or both genes are missing for a FI, that FI will be excluded.
     * @param geneToScore
     * @param useAsPValue the scores are p-values, geometric mean should be used. The edge score
     * should be -log10.
     * @return
     */
    public List<Set<String>> mclClusterForGeneScores(Map<String, Double> geneToScore,
                                                     boolean useAsPValue,
                                                     Integer sizeCutoff,
                                                     Double corCutoff) throws Exception {
        // Load FI network
        Set<String> fis = new FileUtility().loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisWithValues = new HashSet<String>();
        for (String fi : fis) {
            String[] tokens = fi.split("\t");
            Double value1 = geneToScore.get(tokens[0]);
            Double value2 = geneToScore.get(tokens[1]);
            if (value1 == null || value2 == null)
                continue;
            if (useAsPValue) {
                // Use geometric mean for pvalue
                Double score = -Math.log10(Math.sqrt(value1 * value2));
                fisWithValues.add(fi + "\t" + score);
            }
            else {
//                Double score = (value1 + value2) / 2.0d;
                Double score = Math.min(value1, value2);
                fisWithValues.add(fi + "\t" + score);
            }
        }
        
        //        clusterWrapper.setInflation(4.5d);
        List<Set<String>> clusters = mclCluster(fisWithValues);
        filterClusters(clusters, fisWithValues, fis, sizeCutoff, corCutoff);
        return clusters;
    }
    
    public List<Set<String>> mclCluster(Set<String> fisWithCorrs,
                                        Set<String> fis,
                                        int sizeCutoff, 
                                        double corCutoff) throws Exception {
        List<Set<String>> clusters = mclCluster(fisWithCorrs);
        filterClusters(clusters, 
                       fisWithCorrs, 
                       fis, 
                       sizeCutoff,
                       corCutoff);
        return clusters;
    }

    /**
     * Do a filtering based on average score and size.
     * @param clusters
     * @param fisWithCorrs
     * @param fis
     * @param sizeCutoff
     * @param corCutoff
     * @return the same cluster set should be returned.
     */
    private void filterClusters(List<Set<String>> clusters,
                                Set<String> fisWithCorrs,
                                Set<String> fis, 
                                Integer sizeCutoff,
                                Double corCutoff) {
        if (sizeCutoff != null) {
            for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
                Set<String> cluster = it.next();
                if (cluster.size() < sizeCutoff)
                    it.remove();
            }
        }
        if (corCutoff != null) {
            List<Double> averageCorrs = calculateAverageCorrelations(clusters, 
                                                                     fis, 
                                                                     fisWithCorrs);
            // Filter out some clusters
            int index = 0;
            for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
                it.next();
                Double cor = averageCorrs.get(index);
                if (cor == null || cor < corCutoff)
                    it.remove();
                index ++;
            }
        }
    }

    /**
     * The actual method for doing MCL clustering.
     * @param fisWithCorrs
     * @return
     * @throws Exception
     */
    public List<Set<String>> mclCluster(Set<String> fisWithCorrs) throws Exception {
        MCLClusteringResult mclResults = mclHelper.cluster(fisWithCorrs, inflation);
        List<Set<String>> clusters = parseMCLResults(mclResults);
        return clusters;
    }
    
    public List<Double> calculateAverageCorrelations(List<Set<String>> clusters,
                                                     Set<String> fis,
                                                     Set<String> fisWithCorrs) {
        // Use this map for calculation
        Map<String, Double> fiToCorr = new HashMap<String, Double>();
        int index = 0;
        for (String fiWithCorr : fisWithCorrs) {
            index = fiWithCorr.lastIndexOf("\t");
            fiToCorr.put(fiWithCorr.substring(0, index),
                         new Double(fiWithCorr.substring(index + 1)));
        }
        List<Double> corrs = new ArrayList<Double>();
        double total = 0.0d;
        int count = 0;
        index = 0;
        for (Set<String> cluster : clusters) {
            index ++;
            Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
            total = 0.0d;
            count = 0;
            for (String fi : fisInCluster) {
                Double value = fiToCorr.get(fi);
                if (value != null) {
                    total += value;
                    count ++;
                }
            }
            if (count == 0) // It is possible for some small modules or for sparse array data sets
                corrs.add(null);
            else
                corrs.add(total / count);
        }
        return corrs;
    }
    
    private List<Set<String>> parseMCLResults(MCLClusteringResult result) {
        List<Set<String>> rtn = new ArrayList<Set<String>>();
        for (String text : result.getClusters()) {
            String[] tokens = text.split("\t");
            Set<String> set = new HashSet<String>();
            for (String token : tokens)
                set.add(token);
            rtn.add(set);
        }
        return rtn;
    }
    
    
}
