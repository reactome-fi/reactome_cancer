/*
 * Created on Dec 11, 2012
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.inference.TestUtils;
import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do analyses relating exome mutation and gene expression together.
 * @author gwu
 *
 */
public class ExomeMutationAndGeneExpressionCorrelationAnalyzer {
    
    public ExomeMutationAndGeneExpressionCorrelationAnalyzer() {
    }
    
    /**
     * This method is used to get differentially expression genes based on
     * exome mutation file.
     * @throws Exception
     */
    @Test
    public void checkGeneExpressionInDifferentModules() throws Exception {
        Set<String> moduleGenes = loadExomeMutationModuleGenes();
        System.out.println("Module genes: " + moduleGenes.size());
//        for (String gene : moduleGenes)
//            System.out.println(gene);
        
        // The following is just a quick check of the selected module genes.
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        
        MATFileLoader matFileLoader = new MATFileLoader();
        // TCGA OV
//        String dirName = "datasets/TCGA/OvarianCancer/FireHose/";
//        String mutationFileName = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/OV.final_analysis_set.maf";
//        Map<String, Set<String>> sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
//        String clinFileName = dirName + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2012082500.0.0/OV.clin.transformed.txt";
//        String geneExpFile = dirName + "gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_4.2012082500.0.0/OV.medianexp.transformed.global.zscore.txt";
        
        // TCGA GBM
        String dirName = "datasets/TCGA/GBM/FireHose/";
        String mutationFileName = dirName + "FH_GBM.Mutation_Significance.Level_4.20120725/GBM.final_analysis_set.maf";
        Map<String, Set<String>> sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        String clinFileName = dirName + "FH_GBM.Clinical_Pick_Tier1.Level_4.20120825/GBM.clin.transformed.txt";
        String geneExpFile = dirName + "gdac.broadinstitute.org_GBM.mRNA_Preprocess_Median.Level_4.2012102400.0.0/GBM.medianexp.transformed.global.zscore.txt";
        
        Set<String> group1 = new HashSet<String>();
        Set<String> group2 = new HashSet<String>();
        splitSamples(group1, 
                     group2, 
                     sampleToGenes, 
                     moduleGenes);
        System.out.println("Group 1 (no mutation): " + group1.size());
        System.out.println("Group 2 (mutation): " + group2.size());
        TCGAOvarianCancerAnalyzer cancerAnalyzer = new TCGAOvarianCancerAnalyzer();
        cancerAnalyzer.doSurvivalAnalysisForMutationModule(survivalHelper,
                                                           sampleToGenes,
                                                           moduleGenes, 
                                                           clinFileName, 
                                                           "TCGA_OV_Module_");
        
        // Load gene expression data set
        CancerGeneExpressionCommon geneExpressionHelper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = geneExpressionHelper.loadGeneExp(geneExpFile);
        
        // Permutation tests for calculating q-values for the following t-test
        int permutation = 100;
        List<Double> randomPValues = new ArrayList<Double>();
        Set<String> allMutatedGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
        RandomData randomizer = new RandomDataImpl();
        for (int i = 0; i < permutation; i++) {
            System.out.println("Permutation " + i);
            Set<String> randomGenes = MathUtilities.randomSampling(allMutatedGenes,
                                                                   moduleGenes.size(),
                                                                   randomizer);
            Set<String> randomGroup1 = new HashSet<String>();
            Set<String> randomGroup2 = new HashSet<String>();
            splitSamples(randomGroup1, 
                         randomGroup2, 
                         sampleToGenes, 
                         randomGenes);
            for (String gene : geneToSampleToValue.keySet()) {
                Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                double[] ttestResults = doTTest(sampleToValue, 
                                                sampleToGenes.keySet(),
                                                randomGroup1,
                                                gene);
                randomPValues.add(ttestResults[1]);
            }
        }
        Collections.sort(randomPValues);
        
        // Assume all sequenced samples have been listed in this map. This may not be true if some sequenced
        // samples don't have any mutation though this is highly unlikely.
        System.out.println("Total samples in the mutation map: " + sampleToGenes.size());
        // Running a t-test to get a list of top differential expression genes
        List<Double> realPValues = new ArrayList<Double>();
        Map<String, double[]> geneToTTestValues = new HashMap<String, double[]>();
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            double[] ttestResults = doTTest(sampleToValue, 
                                            sampleToGenes.keySet(),
                                            group1,
                                            gene);
            realPValues.add(ttestResults[1]);
            geneToTTestValues.put(gene, ttestResults);
        }
        Collections.sort(realPValues);
        System.out.println("Gene\tT\tP-value\tFDR");
        for (String gene : geneToTTestValues.keySet()) {
            double[] tTestResults = geneToTTestValues.get(gene);
            double fdr = MathUtilities.calculateFDR(tTestResults[1],
                                                    realPValues,
                                                    randomPValues);
            System.out.println(gene + "\t" + tTestResults[0] + "\t" + tTestResults[1] + "\t" + fdr);
        }
    }
    
    private double[] doTTest(Map<String, Double> sampleToValue,
                             Set<String> samples,
                             Set<String> group1,
                             String gene) throws MathException {
        List<Double> list1 = new ArrayList<Double>();
        List<Double> list2 = new ArrayList<Double>();
        for (String sample : samples) {
            Double value = sampleToValue.get(sample);
            if (value == null)
                continue;
//            if (value == null)
//                throw new IllegalStateException(gene + " in " + sample + " has no expression value!");
            if (group1.contains(sample)) // Either in group1 or group2.
                list1.add(value);
            else
                list2.add(value); 
        }
        double[] values1 = convertToArray(list1);
        double[] values2 = convertToArray(list2);
        Double tvalue = TestUtils.t(values1, 
                                    values2);
        Double pvalue = TestUtils.tTest(values1, 
                                        values2);
        return new double[]{tvalue, pvalue};
    }
    
    private double[] convertToArray(List<Double> list) {
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            array[i] = list.get(i);
        return array;
    }
    
    /**
     * Split samples into two groups: one group has genes mutated and another does.
     * @param group1 samples having no genes mutated in moduleGenes
     * @param group2 samples having genes mutated in modulesGenes
     * @param sampleToGenes
     * @param moduleGenes
     */
    private void splitSamples(Set<String> group1,
                              Set<String> group2,
                              Map<String, Set<String>> sampleToGenes,
                              Set<String> moduleGenes) {
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(genes, moduleGenes);
            if (shared.size() > 0)
                group2.add(sample);
            else
                group1.add(sample);
        }
    }
    
    private Set<String> loadExomeMutationModuleGenes() throws IOException {
        // Load TCGA OV mutation genes
        //String moduleFileName = R3Constants.OVARIAN_DIR_NAME + "FireHose/2009FISubNetwork_Sample3_Modules_092012.txt";
        // TCGA GBM mutation genes
        String moduleFileName = R3Constants.GBM_DIR + "FireHose/2012FISubNetwork_Sample5_Modules_092012.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> modules = clusterHelper.loadNetworkClusters(moduleFileName);
        return modules.get(3);
    }
    
}
