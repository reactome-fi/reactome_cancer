/*
 * Created on Apr 29, 2009
 *
 */
package org.reactome.cancer;

import java.awt.Color;
import java.io.IOException;
import java.util.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TestUtils;
import org.gk.util.StringUtils;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.r3.cluster.HierarchicalCluster;
import org.reactome.r3.cluster.HierarchicalClusterNode;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.BreadthFirstSearch.Edge;
import org.reactome.r3.graph.BreadthFirstSearch.TreeNode;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.JGraphTUtilities;
import org.reactome.r3.graph.JungGraphUtilities;
import org.reactome.r3.graph.MCLResultsAnalyzer;
import org.reactome.r3.graph.NetworkModularityCalculator;
import org.reactome.r3.graph.NetworkModule;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import weka.core.Utils;

/**
 * New methods for cancer resequence data analysis.
 * @author wgm
 *
 */
public class NatureGBMAnalyzer extends CancerResequenceDataSetAnalyzer {
    public static final String TCGA_GBM_DIR = "datasets/TCGA/GBM/";
    
    public NatureGBMAnalyzer() {
    }
    
    /**
     * This method is used to compare the network clusters results from two different GBM data sets.
     * @throws IOException
     */
    @Test
    public void compareNetworkClustersFromTwoDatasets() throws Exception {
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String tcgaClusterFile = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        String parsonClusterFile = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        List<Set<String>> tcgaClusters = clusterAnalyzer.loadNetworkClusters(tcgaClusterFile);
        List<Set<String>> parsonsClusters = clusterAnalyzer.loadNetworkClusters(parsonClusterFile);
        // Do a label
        List<String> tcgaClusterLabels = clusterAnalyzer.labelNetworkClusters(tcgaClusters);
        List<String> parsonsClusterLabels = clusterAnalyzer.labelNetworkClusters(parsonsClusters);
        // Do a pair-wise comparison
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("TCGA_Cluster\tLabel\tSize\tParsons_Cluster\tLabel\tSize\tShared\tP-value\tShared_Genes");
        for (int i = 0; i < tcgaClusters.size(); i++) {
            Set<String> tcgaCluster = tcgaClusters.get(i);
            //if (tcgaCluster.size() < 5)
            //    continue;
            String tcgaLabel = tcgaClusterLabels.get(i);
            for (int j = 0; j < parsonsClusters.size(); j++) {
                Set<String> parsonsCluster = parsonsClusters.get(j);
               // if (parsonsCluster.size() < 5)
               //     continue;
                String parsonsLabel = parsonsClusterLabels.get(j);
                Set<String> shared = InteractionUtilities.getShared(tcgaCluster, parsonsCluster);
                if (shared.size() == 0)
                    continue;
                double pvalue = MathUtilities.calculateHypergeometricPValue(fiGenes.size(), 
                                                                            tcgaCluster.size(),
                                                                            parsonsCluster.size(),
                                                                            shared.size());
                System.out.println(i + "\t" + tcgaLabel + "\t" + tcgaCluster.size() + "\t" +
                                   j + "\t" + parsonsLabel + "\t" + parsonsCluster.size() + "\t" +
                                   shared.size() + "\t" + pvalue + "\t" + 
                                   shared);
            }
        }
    }
    
    private Set<String> getSharedGenesInClustersForTwoDatasets() throws IOException {
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String tcgaClusterFile = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        String parsonClusterFile = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        List<Set<String>> tcgaClusters = clusterAnalyzer.loadNetworkClusters(tcgaClusterFile);
        List<Set<String>> parsonsClusters = clusterAnalyzer.loadNetworkClusters(parsonClusterFile);
        // 0 - 0
        Set<String> tcgaCluster = tcgaClusters.get(0);
        Set<String> parsonCluster = parsonsClusters.get(0);
        Set<String> shared = InteractionUtilities.getShared(tcgaCluster, parsonCluster);
        Set<String> totalShared = new HashSet<String>(shared);
        // 1 - 2
        tcgaCluster = tcgaClusters.get(1);
        parsonCluster = parsonsClusters.get(2);
        shared = InteractionUtilities.getShared(tcgaCluster, parsonCluster);
        totalShared.addAll(shared);
        return totalShared;
    }
    
    @Test
    public void checkSharedGenesInClustersFromTwoDatasets() throws Exception {
        Set<String> shared = getSharedGenesInClustersForTwoDatasets();
        System.out.println("Total shared: " + shared.size());
        // Check average shortest path
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<TreeNode, List<Edge>> geneToPartners = bfs.initGraph(fis);
        double averagePath = calculateShortestPath(new ArrayList<String>(shared), 
                                                   bfs, 
                                                   geneToPartners);
        System.out.println("Average dist: " + averagePath);
        // Check how many GBM samples for these genes
        // TCGA data set
        List<String> sequencedSamples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(sequencedSamples);
        final Map<String, Set<String>> geneToTCGASamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        // Parsons data set
        sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        final Map<String, Set<String>> geneToParsonsSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        List<String> geneList = new ArrayList<String>(shared);
        Collections.sort(geneList);
        for (String gene : geneList) {
            Set<String> samples1 = geneToTCGASamples.get(gene);
            Set<String> samples2 = geneToParsonsSamples.get(gene);
            System.out.println(gene + "\t" + samples1.size() + "\t" + samples2.size());
        }
        // This is used to print out the altered genes based on sample numbers
        System.out.println("\nGenes ordered based on samples for TCGA:");
        geneList = new ArrayList<String>(geneToTCGASamples.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Set<String> samples1 = geneToTCGASamples.get(gene1);
                Set<String> samples2 = geneToTCGASamples.get(gene2);
                return samples2.size() - samples1.size();
            }
        });
        for (String gene : geneList) {
            Set<String> samples = geneToTCGASamples.get(gene);
            if (samples.size() < 2)
                continue;
            System.out.println(gene + "\t" + samples.size());
        }
        System.out.println("\nGenes ordered based on samples for Parsons:");
        geneList = new ArrayList<String>(geneToParsonsSamples.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Set<String> samples1 = geneToParsonsSamples.get(gene1);
                Set<String> samples2 = geneToParsonsSamples.get(gene2);
                return samples2.size() - samples1.size();
            }
        });
        for (String gene : geneList) {
            Set<String> samples = geneToParsonsSamples.get(gene);
            if (samples.size() < 2)
                continue;
            System.out.println(gene + "\t" + samples.size());
        }
        //PathwayAnnotator annotator = new PathwayAnnotator();
        //annotator.annotateGenesWithFDR(shared);
    }
    
    
    /**
     * This method is used to check if there is a correlation between the distance and odds ratio
     * between two mutated genes. The odds ratio is used to check condition probability between two 
     * mutated genes to see if there is any conditional relationship between these two genes (e.g. 
     * mutual exclusivity or mutual dependence), which is defined as p(m1, m2) / (p(m1) * p(m2)), p
     * is the probabilty for a mutation in the sample set.
     * @throws IOException
     */
    @Test
    public void checkCorrelationBetweenGeneDistanceAndOddsRatio() throws Exception {
        Map<String, Set<String>> sampleToMutations = loadSampleToMutations();
        int totalSample = sampleToMutations.size();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToMutations);
        // Try a filter
        int sampleCutoff = 1;
        for (Iterator<String> it = geneToSamples.keySet().iterator(); it.hasNext();) {
            Set<String> samples = geneToSamples.get(it.next());
            if (samples.size() < sampleCutoff)
                it.remove();
        }
        System.out.println("Sample cutoff: " + sampleCutoff);
        System.out.println("Total used genes: " + geneToSamples.size());
        // Used to check distance
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        
        List<String> geneList = new ArrayList<String>(geneToSamples.keySet());
        // Want to work genes in our FI network only
        geneList.retainAll(geneToPartners.keySet());
        Collections.sort(geneList);
        List<Double> distanceList = new ArrayList<Double>();
        List<Double> oddsList = new ArrayList<Double>();
        System.out.println("Gene1\tGene2\tDist\tOdds");
        for (int i = 0; i < geneList.size() - 1; i++) {
            String gene1 = geneList.get(i);
            Set<String> samples1 = geneToSamples.get(gene1);
            for (int j = i + 1; j < geneList.size(); j++) {
                String gene2 = geneList.get(j);
                Set<String> samples2 = geneToSamples.get(gene2);
                int dist = bfs.getDistance(gene1, 
                                           gene2, 
                                           geneToPartners);
                distanceList.add((double)dist);
                // Calculate odds
                Set<String> shared = InteractionUtilities.getShared(samples1, samples2);
                double odds = ((double) totalSample) * shared.size() / (samples1.size() * samples2.size());
                oddsList.add(odds);
                System.out.println(gene1 + "\t" + gene2 + "\t" + 
                                   dist + "\t" + odds);
            }
        }
        double cc = MathUtilities.calculatePearsonCorrelation(distanceList, oddsList);
        System.out.println("Pearson correlation: " + cc);
    }
    
    /**
     * This method is used to check mutated genes in MCL clusters.
     * @throws Exception
     */
    @Test
    public void checkAlteredGenesInMCLClusters() throws Exception {
        MCLResultsAnalyzer mclAnalyzer = new MCLResultsAnalyzer();
        List<Set<String>> clusters = mclAnalyzer.loadMCLClusters(R3Constants.RESULT_DIR + "MCLCluster_FIsInGene_041709_I40.txt");
//        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        mutatedGenes.retainAll(allGenes);
        // TCGA GBM data set
        List<String> samples = loadResequencedSamples();
        Set<String> mutatedGenes = getAlteredGenesInSamples(samples);
        mutatedGenes.retainAll(allGenes);
        System.out.println("Total mutated genes: " + mutatedGenes.size());
        int total = allGenes.size();
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            if (cluster.size() < 3)
                continue;
            Set<String> shared = new HashSet<String>(cluster);
            shared.retainAll(mutatedGenes);
            double pvalue = MathUtilities.calculateHypergeometricPValue(total,
                                                                        cluster.size(),
                                                                        mutatedGenes.size(),
                                                                        shared.size());
            System.out.println(i + "\t" +
                               cluster.size() + "\t" +
                               shared.size() + "\t" +
                               pvalue);
        }
    }
    
    @Test
    public void outputAlteredGenes() throws Exception {
        // Output altered genes in a Cytoscape format
        // TCGA
//        List<String> samples = loadResequencedSamples();
//        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        // Parsons
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        //String outFileName = TCGA_GBM_DIR + "TCGAAlteredGenesToSampleNumbers.txt";
        String outFileName = TCGA_GBM_DIR + "ParsonsAlteredGenesToSampleNumbers.txt";
        fu.setOutput(outFileName);
        fu.printLine("Samples (class=java.lang.Integer)");
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples1 = geneToSamples.get(gene);
            fu.printLine(gene + "=" + samples1.size());
        }
        fu.close();
//        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
//        System.out.println("#TCGA GBM Genes with Two or More Somatic Mutations: " + mutatedGenes.size());
//        for (String gene : mutatedGenes)
//            System.out.println(gene);
//        // Altered genes
//        List<String> samples = loadResequencedSamples();
//        Set<String> alteredGenes = getAlteredGenesInSamples(samples);
//        List<String> alteredGeneList = new ArrayList<String>(alteredGenes);
//        System.out.println("\n#TCGA All Somatic and CNV genes: " + alteredGeneList.size());
//        for (String gene : alteredGeneList)
//            System.out.println(gene);    
    }
    
    /**
     * This method is used to check the distribution of samples in the MCL clusters.
     * @throws Exception
     */
    @Test
    public void checkSamplesInMCLClusters() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        //Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        Set<String> allAlteredGenes = getAlteredGenesInSamples(samples);
        MCLResultsAnalyzer mclAnalyzer = new MCLResultsAnalyzer();
        List<Set<String>> clusters = mclAnalyzer.loadMCLClusters(R3Constants.RESULT_DIR + "MCLCluster_FIsInGene_041709_I40.txt");
        System.out.println("Cluster\tCluster_Size\tHit_Samples\tPercentage\tAltered_Gene");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            if (cluster.size() < 3)
                continue;
            //if (i != 33)
            //    continue;
            int hitSamples = 0;
            StringBuilder builder = new StringBuilder();
            for (String sample : sampleToAlteredGenes.keySet()) {
                Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
                Set<String> shared = new HashSet<String>(alteredGenes);
                shared.retainAll(cluster);
                if (shared.size() > 0) {
                    hitSamples ++;
                    //System.out.println(sample + "\t" + shared);
                    builder.append(sample + "\t" + shared).append("\n");
                }
            }
            // Check how many mutated genes
            Set<String> totalShared = new HashSet<String>(cluster);
            totalShared.retainAll(allAlteredGenes);
            //if (hitSamples > 1 && totalShared.size() > 1) {
                System.out.println(i + "\t" + 
                                   cluster.size() + "\t" +
                                   hitSamples + "\t" +
                                   (double)hitSamples / sampleToAlteredGenes.size() + "\t" +
                                   totalShared);
            //    System.out.println(builder.toString());
            //}
        }
    }
    
    /**
     * Use hyper-geometric test to check the top 10 genes in network clusters.
     * @throws Exception
     */
    @Test
    public void checkDistributionOfGenesInClusters() throws Exception {
        // For TCGA genes
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        System.out.println("Top ten genes in TCGA data set:");
        Set<String> tcgaGenes = getAlteredGenesInSamples(samples);
        checkDistributionOfGenesInClusters(sampleToAlteredGenes);
        // For science genes
        System.out.println("\nTop ten genes in Science data set:");
        sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        checkDistributionOfGenesInClusters(sampleToAlteredGenes);
    }
    
    private void checkDistributionOfGenesInClusters(Map<String, Set<String>> sampleToAlteredGenes) throws Exception {
        // Used to filter out genes not in the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(clusterFileName);
        List<String> samples = loadResequencedSamples();
        Set<String> tcgaGenes = getAlteredGenesInSamples(samples);
        // Want to check genes in the FI network only
        tcgaGenes.retainAll(fiGenes);
        int numberInCluster0And1 = clusters.get(0).size() + clusters.get(1).size();
        final Map<String, Set<String>> genesToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        List<String> geneList = new ArrayList<String>(genesToSamples.keySet());
        geneList.retainAll(fiGenes);
        // Want to check tcga genes if the passed genes are from Science GBM.
        geneList.retainAll(tcgaGenes);
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Set<String> samples1 = genesToSamples.get(gene1);
                Set<String> samples2 = genesToSamples.get(gene2);
                return samples2.size() - samples1.size();
            }
        });
        int index = 0;
        for (String gene : geneList) {
            Set<String> set = genesToSamples.get(gene);
            System.out.println(index + "\t" + gene + "\t" + set.size());
            index ++;
        }
        // Check the first ten genes
        int top = 10;
        List<String> topGenes = new ArrayList<String>(geneList.subList(0, top));
        int topInCluster0And1 = 0;
        topGenes.retainAll(clusters.get(0));
        System.out.println("In cluster 0: " + topGenes);
        topInCluster0And1 += topGenes.size();
        topGenes = new ArrayList<String>(geneList.subList(0, top));
        topGenes.retainAll(clusters.get(1));
        System.out.println("In cluster 1: " + topGenes);
        topInCluster0And1 += topGenes.size();
        System.out.println("Total genes in network: " + tcgaGenes.size());
        System.out.println("Network clusters 0 or 1: " + numberInCluster0And1);
        System.out.println("Top genes in clusters 0 or 1: " + topInCluster0And1);
        double pvalue = MathUtilities.calculateHypergeometricPValue(tcgaGenes.size(),
                                                                    top, 
                                                                    numberInCluster0And1, 
                                                                    topInCluster0And1);
        System.out.println("p-value from hyper-geometric: " + pvalue);
    }
    
    /**
     * This method is used to generate a Cytoscape node attribute file for
     * gene to sample numbers.
     * @throws Exception
     */
    @Test
    public void generateGeneToSampleNumberNodeAttFile() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenese = getSampleToAlteredGenes(samples);
        // For TCGA genes
        String fileName = TCGA_GBM_DIR + "TCGAAlteredGenesToSampleNumbers.na";        
        generateGeneToSampleNumberNodeAttFile(fileName, sampleToAlteredGenese);
        // For Science genes
        sampleToAlteredGenese = getScienceGBMSampleToAlteredGenes();
        fileName = TCGA_GBM_DIR + "ScienceAlteredGenesToSampleNumbers.na";
        generateGeneToSampleNumberNodeAttFile(fileName, sampleToAlteredGenese);
    }

    private void generateGeneToSampleNumberNodeAttFile(String fileName,
                                                       Map<String, Set<String>> sampleToAlteredGenese)
            throws IOException {
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenese);
        fu.setOutput(fileName);
        fu.printLine("SampleNumber (class=java.lang.Integer)");
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            fu.printLine(gene + "=" + samples.size());
        }
        fu.close();
    }
    
    @Test
    public void calcualteSampleCorrelationBetweenSurvivalAndDistance() throws Exception {
        List<String> resequencedSamples = loadResequencedSamples();
        //Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(resequencedSamples);
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        String clusterFileName = TCGA_GBM_DIR + "SampleClusteringFromHomoCNVAndMutation_063009.txt";
        // Load sample distance
        Map<String, Double> samplePairToDist = new HashMap<String, Double>();
        fu.setInput(clusterFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            if (line.equals("hierarchical layouting..."))
                break;
            String[] tokens = line.split(": ");
            Double dist = new Double(tokens[1]);
            samplePairToDist.put(tokens[0], dist);
        }
        fu.close();
        Map<String, Integer> sampleToSurvival = loadSampleToTimeSpan();
        double[] distArray = new double[samplePairToDist.size()];
        double[] survivalArray = new double[samplePairToDist.size()];
        double[] sharedGenes = new double[samplePairToDist.size()];
        double[] sharedFIs = new double[samplePairToDist.size()];
        int index = 0;
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Map<String, Set<String>> idToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        for (String samplePair : samplePairToDist.keySet()) {
            Double dist = samplePairToDist.get(samplePair);
            distArray[index] = dist;
            String[] twoSamples = samplePair.split("\t");
            int survival1 = sampleToSurvival.get(twoSamples[0]);
            int survival2 = sampleToSurvival.get(twoSamples[1]);
            survivalArray[index] = Math.abs(survival1 - survival2);
            Set<String> alteredGenes1 = sampleToAlteredGenes.get(twoSamples[0]);
            Set<String> alteredGenes2 = sampleToAlteredGenes.get(twoSamples[1]);
            if (alteredGenes1 == null || alteredGenes2 == null)
                continue;
            double jarccardIndex = MathUtilities.calculateJaccardIndex(alteredGenes1, alteredGenes2); 
            sharedGenes[index] = jarccardIndex;
            Set<String> fis1 = getFIsForGenes(alteredGenes1, idToPartners);
            Set<String> fis2 = getFIsForGenes(alteredGenes2, idToPartners);
            jarccardIndex = MathUtilities.calculateJaccardIndex(fis1, fis2);
            sharedFIs[index] = jarccardIndex;
            index ++;
        }
        double cc = Utils.correlation(distArray, 
                                      survivalArray, 
                                      index);
        System.out.println("Total sample pairs: " + index);
        System.out.println("Correlation coeficient between distances and survival: " + cc);
        cc = Utils.correlation(sharedGenes, survivalArray, survivalArray.length);
        System.out.println("Correlation coeficient between shared genes and survival: " + cc);
        cc = Utils.correlation(sharedFIs, survivalArray, survivalArray.length);
        System.out.println("Correlation coeficient between shared FIs and survival: " + cc);
    }
    
    private Set<String> getFIsForGenes(Set<String> genes,
                                       Map<String, Set<String>> idToPartners) {
        Set<String> rtn = new HashSet<String>();
        int compare = 0;
        for (String gene : genes) {
            Set<String> partners = idToPartners.get(gene);
            if (partners == null)
                continue;
            for (String partner : partners) {
                compare = gene.compareTo(partner);
                if (compare < 0)
                    rtn.add(gene + "\t" + partner);
                else if (compare > 0)
                    rtn.add(partner + "\t" + gene);
            }
        }
        return rtn;
    }
    
    @Test
    public void clusterAlteredGenes() throws Exception {
        // TCGA GBM data set
//        List<String> samples = loadResequencedSamples();
//        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
//        Set<String> alteredGenes = selectGenesInSamples(2, sampleToAlteredGenes);
//        Set<String> alteredGenes = getAlteredGenesInSamples(samples);
        // Science GBM data sets
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        Set<String> alteredGenes = new HashSet<String>();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            alteredGenes.add(gene);
        }
        // Output interactions
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total altered genes: " + alteredGenes.size());
        alteredGenes.retainAll(fiGenes);
        System.out.println("Altered genes in FI network: " + alteredGenes.size());
        Set<String> fisInAlteredGenes = InteractionUtilities.getFIs(alteredGenes, fis);
        Set<String> genesInFIsInAlteredGenes = InteractionUtilities.grepIDsFromInteractions(fisInAlteredGenes);
        System.out.println("Genes connected: " + genesInFIsInAlteredGenes.size());
//        //String outFileName = TCGA_GBM_DIR + "FIsInTCGAAlteredGenes071609_1.txt";
//        String outFileName = TCGA_GBM_DIR + "FIsInScienceGBMAlteredGenes072209.txt";
        //String outFileName = TCGA_GBM_DIR + "FIsInTCGAAlteredGenesInPathwayFIs080509.txt";
//        String outFileName = TCGA_GBM_DIR + "FIsInScienceGBMAlteredGenes092909.txt";
//        fu.saveInteractions(fisInAlteredGenes, outFileName);
//        if (true)
//            return;
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
//        long time1 = System.currentTimeMillis();
        List<Set<String>> clusters = clusterAnalyzer.cluster(alteredGenes, 
                                                             //0.21, // The best ratio
                                                             //0.19, // Best for the Science GBM
                                                             fis);
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for edge betweenness clustering: " + (time2 - time1));
//        //String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
//        //String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes032910.txt";
//        //List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(outFileName);
//        double modularity = clusterAnalyzer.calculateModualarity(clusters, fisInAlteredGenes);
//        System.out.println("Modularity from edge betweenness: " + modularity);
        
        // Spectral partitioning
//        SpectralPartitionNetworkCluster spectralCluster = new SpectralPartitionNetworkCluster();
//        long time11 = System.currentTimeMillis();
//        List<Set<String>> newClusters = spectralCluster.cluster(fisInAlteredGenes);
//        long time21 = System.currentTimeMillis();
//        System.out.println("Time for spectral partitioning: " + (time21 - time11));
//        double modularity = clusterAnalyzer.calculateModularity(newClusters, 
//                                                          fisInAlteredGenes);
//        System.out.println("Modularity from spectral partition: " + modularity);
//        String outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes042810.txt";
        String outFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes042810.txt";
        clusterAnalyzer.outputNetworkClusters(clusters, outFileName);
        
        int index = 0;
        Set<String> totalClusterGenes = new HashSet<String>();
        for (Set<String> cluster : clusters) {
            int size = cluster.size();
            double percent = (double) size / alteredGenes.size();
            System.out.println(index + "\t" +  + size + "\t" + percent);
            index++;
            totalClusterGenes.addAll(cluster);
        }
//        System.out.println("Total genes in cluster: " + totalClusterGenes.size());
//        String outFileName = TCGA_GBM_DIR + "ClustersInParsonsAlteredGenes_CandidateCNV_093009.txt";
//        outFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesInPathwayFIs080509_220.txt";
//        clusterAnalyzer.outputNetworkClusters(clusters, outFileName);
    }
    
    @Test
    public void compareTwoNetworkClusters() throws IOException {
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String fileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        //String fileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes042810.txt";
        List<Set<String>> scienceClusters = clusterAnalyzer.loadNetworkClusters(fileName);
        //fileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        fileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes042810.txt";
        List<Set<String>> tcgaClusters = clusterAnalyzer.loadNetworkClusters(fileName);
        // Check sharing
        System.out.println("TCGA_Cluster\tSize\tScience_Cluster\tSize\tShared\tShared_Percentage");
        // Just want to see the first five clusters
        for (int i = 0; i < tcgaClusters.size(); i++) {
            Set<String> cluster1 = tcgaClusters.get(i);
            for (int j = 0; j < scienceClusters.size(); j++) {
                Set<String> cluster2 = scienceClusters.get(j);
                //System.out.println(j + ": " + cluster2);
                int min = Math.min(cluster1.size(), cluster2.size());
                Set<String> copy = new HashSet<String>(cluster2);
                copy.retainAll(cluster1);
                double shared = copy.size() / (double) min;
                if (copy.size() == 0)
                    continue;
                System.out.println(i + "\t" + cluster1.size() + "\t" +
                                   j + "\t" + cluster2.size() + "\t" +
                                   copy.size() + "\t" + shared);
            }
        }
    }
    
    @Test
    public void clusterTwoGBMAlteredGenes() throws Exception {
        Map<String, Set<String>> tcgaSampleToMutatedGenes = loadSampleToMutations();
        Map<String, Set<String>> scienceSampleToMutatedGenes = loadAllScienceSampleToMutatedGenes();
        Set<String> tcgaMutatedGenesIn2Samples = selectGenesInSamples(2, tcgaSampleToMutatedGenes);
        Set<String> scienceMutatedGenesIn2Samples = selectGenesInSamples(2, scienceSampleToMutatedGenes);
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> tcgaClusters = clusterAnalyzer.cluster(tcgaMutatedGenesIn2Samples);
        String fileName = TCGA_GBM_DIR + "TCGANetworkClustersForMutatedGenesIn2Samples.txt";
        clusterAnalyzer.outputNetworkClusters(tcgaClusters, fileName);
        List<Set<String>> scienceClusters = clusterAnalyzer.cluster(scienceMutatedGenesIn2Samples);
        fileName = TCGA_GBM_DIR + "ScienceNetworkClustersForMutatedGenesIn2Samples.txt";
        clusterAnalyzer.outputNetworkClusters(scienceClusters, fileName);
        // Want to print FIs
        Set<String> allFIs = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> tcgaFIs = InteractionUtilities.getFIs(tcgaMutatedGenesIn2Samples, 
                                                          allFIs);
        fileName = TCGA_GBM_DIR + "FIsInTCGAMutatedGenesIn2Samples.txt";
        fu.saveInteractions(tcgaFIs, fileName);
        Set<String> scienceFIs = InteractionUtilities.getFIs(scienceMutatedGenesIn2Samples,
                                                             allFIs);
        fileName = TCGA_GBM_DIR + "FIsInScienceMutatedGenesIn2Samples.txt";
        fu.saveInteractions(scienceFIs, fileName);
        // Try to merge these two gene sets together
        Set<String> merged = new HashSet<String>();
        merged.addAll(tcgaMutatedGenesIn2Samples);
        merged.addAll(scienceMutatedGenesIn2Samples);
        Set<String> mergedFIs = InteractionUtilities.getFIs(merged, allFIs);
        fileName = TCGA_GBM_DIR + "FIsInTwoGBMMutatedGenesIn2Samples.txt";
        fu.saveInteractions(mergedFIs, fileName);
        // Cluster for the merged FIs
        List<Set<String>> mergedClusters = clusterAnalyzer.cluster(merged);
        fileName = TCGA_GBM_DIR + "TwoGBMNetworkClustersForMutatedGenesIn2Samplex.txt";
        clusterAnalyzer.outputNetworkClusters(mergedClusters, fileName);
        System.out.println("\nCalculating modularity:");
        double modularity = clusterAnalyzer.calculateModularity(mergedClusters, mergedFIs);
        System.out.println("Merged: " + modularity);
        modularity = clusterAnalyzer.calculateModularity(tcgaClusters, tcgaFIs);
        System.out.println("TCGA clusters: " + modularity);
        modularity = clusterAnalyzer.calculateModularity(scienceClusters, scienceFIs);
        System.out.println("Science clusters: " + modularity);
    }
    
    /**
     * This method is used to analyze that clustering can enrich the numbers of shared genes.
     * @throws Exception
     */
    @Test
    public void checkOverlappingOfTwoGBMClusters() throws Exception {
        String clusterFileName = TCGA_GBM_DIR + "ClusterListForParsonsGBM100909.txt";
        String cutoff = "3.2170542635658914";
        List<String> scienceCluster = getClusterForAnalysis(clusterFileName, cutoff);
        clusterFileName = TCGA_GBM_DIR + "ClusterListForTCGAGBM100909.txt";
        cutoff = "2.5245098039215685";
        List<String> natureCluster = getClusterForAnalysis(clusterFileName, cutoff);
        System.out.println("Check overlapping in the two clusters:");
        checkOverlappingOfTwoGBMDatasets(scienceCluster, 
                                         natureCluster);
        // Check the original shared genes
        //List<String> scienceGenes = loadScienceGBMGenes();
        //List<String> natureGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        Set<String> geneSet = selectGenesInSamples(2, sampleToAlteredGenes);
        List<String> scienceGenes = new ArrayList<String>(geneSet);
        
        List<String> samples = loadResequencedSamples();
        sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        geneSet = selectGenesInSamples(2, sampleToAlteredGenes);
        List<String> natureGenes = new ArrayList<String>(geneSet);
        // Check how many genes are shared in these two sets
        List<String> copy = new ArrayList<String>(scienceGenes);
        copy.retainAll(natureGenes);
        System.out.println("\nShared all: " + copy.size());
        System.out.println("\nCheck overlapping for all candidate genes:");
        checkOverlappingOfTwoGBMDatasets(scienceGenes, natureGenes);
        // Did a permutation test for the enrichment by clustering
        int target = 12;
        int permutation = 10000;
        int count = 0;
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        scienceGenes.retainAll(fiGenes);
        natureGenes.retainAll(fiGenes);
        RandomData randomizer = new RandomDataImpl();
        // Test random sampling
//        for (int i = 0; i < 10; i++) {
//            Set<String> sample = MathUtilities.randomSampling(scienceGenes,
//                                                              scienceCluster.size(),
//                                                              randomizer);
//            List<String> list = new ArrayList<String>(sample);
//            Collections.sort(list);
//            System.out.println(list.size() + ": " + list);
//        }
        for (int i = 0; i < permutation; i++) {
            Set<String> scienceSample = MathUtilities.randomSampling(scienceGenes, 
                                                                     scienceCluster.size(),
                                                                     randomizer);
            //System.out.println("Science sample: " + scienceSample);
            Set<String> natureSample = MathUtilities.randomSampling(natureGenes, 
                                                                    natureCluster.size(),
                                                                    randomizer);
            //System.out.println("Nature sample: " + natureSample);
            // Check how many genes are shared
            scienceSample.retainAll(natureSample);
            if (scienceSample.size() >= target)
                count ++;
        }
        System.out.println("\nPermutation test result:");
        System.out.println("pvalue: " + (double) count / permutation);
    }
    
    @Test
    public void checkOverlappingOfAlteredGenesInTwoGBMDatasets() throws Exception {
        // For TCGA genes
        List<String> resequencedSamples = loadResequencedSamples();
        Map<String, Set<String>> tcgaSampleToAlteredGenes = getSampleToAlteredGenes(resequencedSamples);
        Set<String> tcgaAlteredGenes = getAlteredGenesInSamples(resequencedSamples);
        // For Science genes
        Map<String, Set<String>> scienceSampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        Set<String> scienceAlteredGenes = new HashSet<String>();
        for (Set<String> set : scienceSampleToAlteredGenes.values())
            scienceAlteredGenes.addAll(set);
        // Check overlapping
        System.out.println("Overlapping of altered genes:");
        checkOverlappingOfTwoGBMDatasets(tcgaAlteredGenes, 
                                         scienceAlteredGenes);
        // Want to check the distribution of shared genes in cluster 0 and 1
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        tcgaAlteredGenes.retainAll(fiGenes);
        int totalTCGAInFIs = tcgaAlteredGenes.size();
        Set<String> shared = new HashSet<String>(tcgaAlteredGenes);
        shared.retainAll(scienceAlteredGenes);
        int sharedNumber = shared.size();
        // Limit the shared genes into the cluster 0, 1
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String fileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(fileName);
        Set<String> twoClusterGenes = new HashSet<String>();
        twoClusterGenes.addAll(clusters.get(0));
        twoClusterGenes.addAll(clusters.get(1));
        tcgaAlteredGenes.retainAll(twoClusterGenes);
        scienceAlteredGenes.retainAll(twoClusterGenes);
        System.out.println("\nOverlapping genes in the first two clusters:");
        shared.retainAll(twoClusterGenes);
        System.out.println("Total TCGA altered genes in FI: " + totalTCGAInFIs);
        System.out.println("Two clustered genes: " + twoClusterGenes.size());
        System.out.println("Shared genes in FI: " + sharedNumber);
        System.out.println("Shared genes in two clusters: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalTCGAInFIs,
                                                                    sharedNumber,
                                                                    twoClusterGenes.size(),
                                                                    shared.size());
        System.out.println("pvalue: " + pvalue);
    }

    @Test
    public void checkOverlappingOfMutatedGenesInTwoGBMDatasets2() throws Exception {
        Map<String, Set<String>> tcgaSampleToGenes = loadSampleToMutations();
        // Get mutated genes
        Set<String> tcgaGenes = new HashSet<String>();
        for (String sample : tcgaSampleToGenes.keySet()) {
            Set<String> genes = tcgaSampleToGenes.get(sample);
            tcgaGenes.addAll(genes);
        }
        Set<String> scienceGenes = loadAllScienceGBMMutatedGenes();
        
        // Check more than 2 sample genes
        Map<String, Set<String>> scienceSampleToGenes = loadAllScienceSampleToMutatedGenes();
        int sampleCutoff = 2;
        Set<String> tcgaGenesIn2Samples = selectGenesInSamples(2, tcgaSampleToGenes);
        System.out.println("\nTwo sample mutated genes in TCGA: " + tcgaGenesIn2Samples.size());
        Set<String> scienceGenesIn2Samples = selectGenesInSamples(2, scienceSampleToGenes);
        System.out.println("Two sample mutated genes in Science: " + scienceGenesIn2Samples.size());
        checkOverlappingOfTwoGBMDatasets(tcgaGenesIn2Samples,
                                         scienceGenesIn2Samples);
        // Try to merge two data sets together as one.
        Map<String, Set<String>> mergedSampleToGenes = new HashMap<String, Set<String>>();
        mergedSampleToGenes.putAll(scienceSampleToGenes);
        mergedSampleToGenes.putAll(tcgaSampleToGenes);
        System.out.println("\nCreating a merged dat set:");
        System.out.println("TCGA Samples: " + tcgaSampleToGenes.size());
        System.out.println("Science Samples: " + scienceSampleToGenes.size());
        System.out.println("Total samples: " + mergedSampleToGenes.size());
        Set<String> mergedGenes = selectGenesInSamples(2, mergedSampleToGenes);
        System.out.println("Two sample mutated genes in TCGA and Science merged: " + mergedGenes.size());
        // Check if FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        mergedGenes.retainAll(fiGenes);
        System.out.println("In FI network: " + mergedGenes.size());
        Set<String> fisInMerged = InteractionUtilities.getFIs(mergedGenes, fis);
        String outFileName = TCGA_GBM_DIR + "FIsInMutatedGenesFromTwoGBMs071509.txt";
        String clusterOutFileName = TCGA_GBM_DIR + "ClustersFIsInMutatedGenesFromTwoGBMs071509_010.txt";
        fu.saveInteractions(fisInMerged, outFileName);
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterAnalyzer.cluster(mergedGenes, 0.10);
        Set<String> genesInClusters = new HashSet<String>();
        for (Set<String> cluster : clusters)
            genesInClusters.addAll(cluster);
        System.out.println("Genes in clusters: " + genesInClusters.size());
        clusterAnalyzer.outputNetworkClusters(clusters, clusterOutFileName);
        // Check genes from the TCGA in the merged set
//        Set<String> tcgaGenesCopy = new HashSet<String>(tcgaGenes);
//        tcgaGenesCopy.retainAll(mergedGenes);
//        // Check genes from the Science in the merged set
//        Set<String> scienceGenesCopy = new HashSet<String>(scienceGenes);
//        scienceGenesCopy.retainAll(mergedGenes);
//        // Check shared in these two copies
//        System.out.println("\nCheck overlapping for merged genes:");
//        checkOverlappingOfTwoGBMDatasets(tcgaGenesCopy, scienceGenesCopy);
//        Set<String> shared = new HashSet<String>(tcgaGenesCopy);
//        shared.retainAll(scienceGenesCopy);
//        // Create a node attribute
//        String gbmSourceFile = TCGA_GBM_DIR + "SourcesForFIsInMutatedGenesFromTwoGBMs071509.txt";
//        fu.setOutput(gbmSourceFile);
//        fu.printLine("GBM_Source (class=java.lang.String)");
//        for (String gene : mergedGenes) {
//            if (shared.contains(gene)) {
//                fu.printLine(gene + "=TCGA,Science");
//            }
//            else if (tcgaGenesCopy.contains(gene))
//                fu.printLine(gene + "=TCGA");
//            else if (scienceGenesCopy.contains(gene))
//                fu.printLine(gene + "=Science");
//        }
//        fu.close();
//        // Check FIs
//        Set<String> tcgaFIsInMerged = InteractionUtilities.getFIs(tcgaGenesCopy, fis);
//        Set<String> genesInTcgaFIsInMerged = InteractionUtilities.grepIDsFromInteractions(tcgaFIsInMerged);
//        System.out.println("TCGA genes in FIs in merged: " + genesInTcgaFIsInMerged.size());
//        Set<String> scienceFIsInMerged = InteractionUtilities.getFIs(scienceGenesCopy, fis);
//        Set<String> genesInScienceFIsInMerged = InteractionUtilities.grepIDsFromInteractions(scienceFIsInMerged);
//        System.out.println("Science genes in FIs in merged: " + genesInScienceFIsInMerged.size());
    }
    
    /**
     * Analyze the mutated genes only.
     * @throws Exception
     */
    @Test
    public void checkOverlappingOfMutatedGenesInTwoGBMDataSets() throws Exception {
        //Map<String, Set<String>> tcgaSampleToGenes = loadSampleToMutations();
        List<String> resequencedSamples = loadResequencedSamples();
        Map<String, Set<String>> tcgaSampleToGenes = getSampleToAlteredGenes(resequencedSamples);
        // Get mutated genes
        Set<String> tcgaGenes = new HashSet<String>();
        for (String sample : tcgaSampleToGenes.keySet()) {
            Set<String> genes = tcgaSampleToGenes.get(sample);
            tcgaGenes.addAll(genes);
        }
        Set<String> scienceGenes = loadAllScienceGBMMutatedGenes();
        checkOverlappingOfTwoGBMDatasets(tcgaGenes, scienceGenes);
        // Check more than 2 sample genes
        //Map<String, Set<String>> scienceSampleToGenes = loadAllScienceSampleToMutatedGenes();
        Map<String, Set<String>> scienceSampleToGenes = getScienceGBMSampleToAlteredGenes();
        int sampleCutoff = 2;
        Set<String> tcgaGenesIn2Samples = selectGenesInSamples(2, tcgaSampleToGenes);
        System.out.println("\nTwo sample mutated genes in TCGA: " + tcgaGenesIn2Samples);
        Set<String> scienceGenesIn2Samples = selectGenesInSamples(2, scienceSampleToGenes);
        System.out.println("Two sample mutated genes in Science: " + scienceGenesIn2Samples);
        // Get not shared genes
        Set<String> notShared = new HashSet<String>(scienceGenesIn2Samples);
        notShared.removeAll(tcgaGenesIn2Samples);
        System.out.println("Not shared: " + notShared.size() + ": " + notShared);
        System.out.println("\nIn two or more samples:");
        checkOverlappingOfTwoGBMDatasets(tcgaGenesIn2Samples, 
                                         scienceGenesIn2Samples);
        // Want to check FI partners
        System.out.println("\nChecking for FI partners in the Science mutated gene sets:");
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        System.out.println("\nCheck all mutated genes:");
        checkFIOverlappingOfTwoGBMDatasets(tcgaGenes, 
                                           scienceGenes, 
                                           geneToPartners);
        System.out.println("\nCheck mutated genes in two or more samples:");
        checkFIOverlappingOfTwoGBMDatasets(tcgaGenesIn2Samples, 
                                           scienceGenesIn2Samples,
                                           geneToPartners);
    }
    
    private void checkFIOverlappingOfTwoGBMDatasets(Set<String> tcgaGenes,
                                                    Set<String> scienceGenes,
                                                    Map<String, Set<String>> geneToPartners) throws MathException {
        Set<String> tcgaFIPartners = new HashSet<String>();
        for (String gene : tcgaGenes) {
            tcgaFIPartners.addAll(geneToPartners.get(gene));
        }
        tcgaFIPartners.removeAll(tcgaGenes);
        System.out.println("Total TCGA FI partners: " + tcgaFIPartners.size());
        Set<String> shared = new HashSet<String>(scienceGenes);
        shared.retainAll(tcgaFIPartners);
        System.out.println("FI partners in Science genes: " + shared.size());
        Set<String> totalGenes = new HashSet<String>(geneToPartners.keySet());
        totalGenes.removeAll(tcgaGenes);
        Set<String> scienceGenesCopy = new HashSet<String>(scienceGenes);
        scienceGenesCopy.removeAll(tcgaGenes);
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(), 
                                                                    scienceGenesCopy.size(),
                                                                    tcgaFIPartners.size(), 
                                                                    shared.size());
        System.out.println("pvalue from hypergeometric test: " + pvalue);
    }

    private void checkOverlappingOfTwoGBMDatasets(Collection<String> tcgaGenes,
                                                  Collection<String> scienceGenes) throws IOException, MathException {
        System.out.println("Total tcga mutated genes: " + tcgaGenes.size());
        System.out.println("Total science GBM genes: " + scienceGenes.size());
        // Check shared genes
        Set<String> sharedGenes = new HashSet<String>(scienceGenes);
        sharedGenes.retainAll(tcgaGenes);
        System.out.println("Total shared: " + sharedGenes.size());
        // Check these genes in the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        tcgaGenes.retainAll(fiGenes);
        scienceGenes.retainAll(fiGenes);
        sharedGenes.retainAll(fiGenes);
        System.out.println("Genes in FI network: " + fiGenes.size());
        System.out.println("TCGA genes: " + tcgaGenes.size());
        System.out.println("Science genes: " + scienceGenes.size());
        System.out.println("Shared: " + sharedGenes.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(fiGenes.size(),
                                                                    tcgaGenes.size(),
                                                                    scienceGenes.size(), 
                                                                    sharedGenes.size());
        System.out.println("pvalue from hyper-geometic test: " + pvalue);
        System.out.println("Shared genes: " + sharedGenes);
    }
    
    /**
     * This method is used to check the cluster from hierarchical clustering onto the clusters from
     * EdgeBetweenness clusters to see if the cluster from the first method is split in the first two 
     * clusters from edge betweenness clustering.
     * @throws IOException
     */
    @Test
    public void checkHierarchicalClusterInEdgeBetweennessClusters() throws IOException {
        String clusterFileName = DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt";
        String cutoff = "1.0";
        List<String> hierarhicalClusterGenes = getClusterForAnalysis(clusterFileName, cutoff);
        System.out.println("Total genes from hierarchical clustering: " + hierarhicalClusterGenes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        List<Set<String>> edgeBetweennessClusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        System.out.println("index\tsize\thierarchical_genes");
        for (int i = 0; i < edgeBetweennessClusters.size(); i++) {
            Set<String> cluster = edgeBetweennessClusters.get(i);
            Set<String> copy = new HashSet<String>(hierarhicalClusterGenes);
            copy.retainAll(cluster);
            System.out.println(i + "\t" + cluster.size() + "\t" + copy.size());
        }
    }
    
    /**
     * This method is used to check the distribution of altered genes in 
     * pathway FIs and predicted FIs.
     * @throws Exception
     */
    @Test
    public void checkAlteredGenesInPathways() throws Exception {
        List<String> samples = loadResequencedSamples();
        Set<String> alteredGenes = getAlteredGenesInSamples(samples);
        // Check these genes in pathways
        // Check these genes in FIs
        Set<String> pathwayFIs = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        Set<String> predictedFIs = fu.loadInteractions(R3Constants.RESULT_DIR + "FIsInGene_Predicted_041709.txt");
        Set<String> predictedGenes = InteractionUtilities.grepIDsFromInteractions(predictedFIs);
        predictedGenes.removeAll(pathwayGenes);
        System.out.println("Pathway genes: " + pathwayGenes.size());
        System.out.println("Predicted genes: " + predictedGenes.size());
        // Check the distribution of altered genes in these two samples
        Set<String> copy = new HashSet<String>(alteredGenes);
        copy.retainAll(pathwayGenes);
        System.out.println("Altered genes in pathways: " + copy.size());
        int size1 = copy.size();
        copy = new HashSet<String>(alteredGenes);
        copy.retainAll(predictedGenes);
        System.out.println("Altered genes in predicted: " + copy.size());
        int size2 = copy.size();
        double zvalue = MathUtilities.calculateZValue(pathwayGenes.size(), pathwayGenes.size() + predictedGenes.size(), 
                                                      size1, size1 + size2);
        System.out.println("Z value: " + zvalue);
        double pvalue = MathUtilities.calTwoTailStandardNormalPvalue(zvalue);
        System.out.println("Two tailed pvalue: " + pvalue);
        double percentage = (double) size1 / (size1 + size2);
        System.out.println("percentage: " + percentage);
    }
    
    @Test
    public void outputHierarchicalClusters() throws IOException {
        List<String> samples = loadResequencedSamples();
        String fileName = TCGA_GBM_DIR + "OneSampleHierarchicalClusterFromNetworkClusters_1.txt";
        Set<String> clusters = fu.loadInteractions(fileName);
        int index = 1;
        for (String cluster : clusters) {
            String[] tokens = cluster.split(", ");
            for (String sample : tokens) {
                System.out.println(sample + "\t" + index);
                samples.remove(sample);
            }
            index ++;
        }
        for (String sample : samples)
            System.out.println(sample + "\t" + index);
    }
    
    public Map<String, List<Integer>> loadSampleToNetworkClusters() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        String networkClusterName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        List<Set<String>> networkClusters = new NetworkClusterAnalyzer().loadNetworkClusters(networkClusterName);
        Map<String, List<Integer>> sampleToNetworkClusters = new HashMap<String, List<Integer>>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            List<Integer> clusters = new ArrayList<Integer>();
            for (int i = 0; i < networkClusters.size(); i++) {
                Set<String> cluster = networkClusters.get(i);
                Set<String> shared = InteractionUtilities.getShared(alteredGenes, cluster);
                if (shared.size() > 0)
                    clusters.add(i);
            }
            sampleToNetworkClusters.put(sample, clusters);
        }
        return sampleToNetworkClusters;
    }
    
    @Test
    public void checkSampleRecurrentInHierarchicalClusters() throws Exception {
        // Check features in sample clusters
        //String clusterFileName = TCGA_GBM_DIR + "SampleClusteringFromHomoCNVAndMutation_063009.txt";
        String clusterFileName = TCGA_GBM_DIR + "TCGASampleClustersFromNetworkClusters_Weighted.txt";
        List<HierarchicalClusterNode> sampleClusters = new HierarchicalCluster().loadHierarchicalClusters(clusterFileName);
        Map<String, String> sampleToRecurrent = loadSampleToRecurrent();
        // Get the total
        int totalNo = 0;
        int totalRec = 0;
        for (String sample : sampleToRecurrent.keySet()) {
            String recurrent = sampleToRecurrent.get(sample);
            if (recurrent.equals("No"))
                totalNo ++;
            else
               totalRec ++;
        }
        FisherExact fisher = new FisherExact(200);
        for (HierarchicalClusterNode cluster : sampleClusters) {
            if (cluster.ids.size() < 2)
                continue;
            int clusterNo = 0;
            int clusterRec = 0;
            for (String sample : cluster.ids) {
                String recurrent = sampleToRecurrent.get(sample);
                if (recurrent.equals("No"))
                    clusterNo ++;
                else
                    clusterRec ++;
            }
            double pvalue = fisher.getTwoTailedP(clusterNo, 
                                                 clusterRec,
                                                 (totalNo - clusterNo),
                                                 (totalRec - clusterRec));
            if (pvalue > 0.10)
                continue;
            System.out.println(cluster.pathDistance + " : " + cluster.ids.size() + " : " + cluster.ids);
            if (cluster.ids.size() == 58) {
                for (String id : cluster.ids)
                    System.out.println(id + "\t1");
                List<String> sampleList = new ArrayList<String>(sampleToRecurrent.keySet());
                sampleList.removeAll(cluster.ids);
                for (String id : sampleList)
                    System.out.println(id + "\t2");
            }
            // Matrix for Fisher test
            System.out.println("Matrix: " + clusterNo + ", " + clusterRec + ", " + 
                               (totalNo - clusterNo) + ", " + (totalRec - clusterRec));
            System.out.println("pvalue from Fisher: " + pvalue);
            System.out.println();
        }
    }

    
    @Test
    public void checkSampleSurvivalInHierarchicalClusters() throws Exception {
        Map<String, String> sampleToRecurrence = loadSampleToRecurrent();
        List<String> primarySamples = new ArrayList<String>();
        for (String sample : sampleToRecurrence.keySet()) {
            String recurrence = sampleToRecurrence.get(sample);
            if (recurrence.equals("No")) 
                primarySamples.add(sample);
        }
        // Check features in sample clusters
        String clusterFileName = TCGA_GBM_DIR + "SampleClusteringFromHomoCNVAndMutation_063009.txt";
        //String clusterFileName = TCGA_GBM_DIR + "TCGASampleClustersFromNetworkClusters.txt";
        List<HierarchicalClusterNode> sampleClusters = new HierarchicalCluster().loadHierarchicalClusters(clusterFileName);
        final Map<String, Integer> sampleToSurvival = loadSampleToSurvivalRate();
        sampleToSurvival.keySet().retainAll(primarySamples);
        // Get the average
        double total = 0.0;
        for (String sample : sampleToSurvival.keySet()) {
            String recurrence = sampleToRecurrence.get(sample);
            Integer rate = sampleToSurvival.get(sample);
            total += rate;
        }
        double mean = total / sampleToSurvival.size();
        System.out.println("Mean: " + mean);
        int index = 0;
        for (HierarchicalClusterNode cluster : sampleClusters) {
            cluster.ids.retainAll(primarySamples);
            if (cluster.ids.size() < 2)
                continue;
            // Want to compare ttest
            double[] values = new double[cluster.ids.size()];
            int i = 0;
            total = 0.0;
            for (String sample : cluster.ids) {
                String tmp = sample.substring(0, 12);
                Integer span = sampleToSurvival.get(tmp);
                values[i] = (double) span;
                i ++;
            }
            double pvalue = TestUtils.tTest(mean, values);
            if (pvalue > 0.50)
                continue;
            //if (cluster.pathDistance > 0.0 || cluster.ids.size() < 5)
            //    continue;
            //if (pvalue > 0.1 || cluster.ids.size() < 5)
            //    continue;
            //if (cluster.pathDistance > 1.0)
            //    continue;
            // calculate the distribution
            System.out.println(cluster.pathDistance + ": " + cluster.ids.size() + ": " + cluster.ids);
            System.out.println("P value from ttest: " + pvalue);
            System.out.println();
            index ++;
        }
    }

    private void checkSampleSurvivalInSampleClusters(List<HierarchicalClusterNode> sampleClusters)
            throws IOException, MathException {
        final Map<String, Integer> sampleToSurvival = loadSampleToSurvivalRate();
        // Get the average
        double total = 0.0;
        for (String sample : sampleToSurvival.keySet()) {
            Integer rate = sampleToSurvival.get(sample);
            total += rate;
        }
        double mean = total / sampleToSurvival.size();
        System.out.println("Mean: " + mean);
        int index = 0;
        for (HierarchicalClusterNode cluster : sampleClusters) {
            if (cluster.ids.size() < 2)
                continue;
            // Want to compare ttest
            double[] values = new double[cluster.ids.size()];
            int i = 0;
            total = 0.0;
            for (String sample : cluster.ids) {
                String tmp = sample.substring(0, 12);
                Integer span = sampleToSurvival.get(tmp);
                values[i] = (double) span;
                i ++;
            }
            double pvalue = TestUtils.tTest(mean, values);
            //if (cluster.pathDistance > 0.0 || cluster.ids.size() < 5)
            //    continue;
            //if (pvalue > 0.1 || cluster.ids.size() < 5)
            //    continue;
            // calculate the distribution
            System.out.println(cluster.pathDistance + ": " + cluster.ids.size() + ": " + cluster.ids);
            System.out.println("P value from ttest: " + pvalue);
            index ++;
        }
    }
    
    @Test
    public void checkSampleFeatures() throws Exception {
        final Map<String, Integer> sampleToTimeSpan = loadSampleToTimeSpan();
        List<String> samples = new ArrayList<String>(sampleToTimeSpan.keySet());
        Collections.sort(samples, new Comparator<String>() {
            public int compare(String sample1, String sample2) {
                Integer span1 = sampleToTimeSpan.get(sample1);
                Integer span2 = sampleToTimeSpan.get(sample2);
                return span2.compareTo(span1);
            }
        });
        List<String> checkingSamples = loadResequencedSamples();
        int c = 1;
        // Create class bins for samples
        int[] groups = new int[2];
        int cutoff = 2;
        for (String sample : samples) {
            if (checkingSamples.contains(sample)) {
                //System.out.println(c + "\t" + sample + "\t" + sampleToTimeSpan.get(sample));
                c ++;
                int span = sampleToTimeSpan.get(sample);
                if (span < cutoff)
                    groups[0] ++;
                else
                    groups[1] ++;
            }
        }
        System.out.println("\tshort\tlong");
        System.out.println("All\t" + groups[0] + "\t" + groups[1]);
        // Check features in sample clusters
        String clusterFileName = TCGA_GBM_DIR + "SampleClusteringFromHomoCNVAndMutation_063009.txt";
        //String clusterFileName = TCGA_GBM_DIR + "SampleClustersUsingHierchicalNetworkCluster.txt";
        //String clusterFileName = TCGA_GBM_DIR + "SampleClustersAvgShortestPath070809.txt";
        //List<ClusterNode> sampleClusters = new HierarchicalClustering().loadHierarchicalClusters(clusterFileName, 1.75);
        List<HierarchicalClusterNode> sampleClusters = new HierarchicalCluster().loadHierarchicalClusters(clusterFileName);
        int index = 0;
        for (HierarchicalClusterNode cluster : sampleClusters) {
            if (cluster.ids.size() < 2)
                continue;
            // calculate the distribution
            System.out.println(cluster.pathDistance + ": " + cluster.ids);
            int[] subgroups = new int[2];
            for (String sample : cluster.ids) {
                String tmp = sample.substring(0, 12);
                Integer span = sampleToTimeSpan.get(tmp);
                if (span == null) {
                    //System.out.println(tmp + " has no time span!");
                    continue;
                }
                if (span < cutoff)
                    subgroups[0] ++;
                else
                    subgroups[1] ++;
            }
            if (subgroups[0] == 0 &&
                subgroups[1] == 0)
                continue;
            System.out.println(index + "\t" + subgroups[0] + "\t" + subgroups[1]);
            // Hypergeometric test should be used since this is no-replacement picking
            double pvalue = MathUtilities.calculateHypergeometricPValue(groups[0] + groups[1],
                                                                        groups[0],
                                                                        subgroups[0] + subgroups[1],
                                                                        subgroups[0]);
            System.out.println("pvalue from hypergeom test for enrichment of the first category: " + pvalue);
            pvalue = MathUtilities.calculateHypergeometricPValue(groups[0] + groups[1],
                                                                 groups[1],
                                                                 subgroups[0] + subgroups[1],
                                                                 subgroups[1]);
            System.out.println("pvalue from hypergeom test for enrichment of the second category: " + pvalue);
            index ++;
        }
    }

    public Map<String, Integer> loadSampleToTimeSpan() throws IOException {
        String fileName = TCGA_GBM_DIR + "IndividualSamples.txt";
        fu.setInput(fileName);
        String line = null;
        // Escape three lines
        for (int i = 0; i < 3; i++)
            line = fu.readLine();
        final Map<String, Integer> sampleToTimeSpan = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[0];
            String start = tokens[6];
            if (start.equals("UNK"))
                continue;
            String end = tokens[7];
            if (end.equals("UNK"))
                continue;
            Integer span = Integer.parseInt(end) - Integer.parseInt(start);
            sampleToTimeSpan.put(sample, span);
        }
        fu.close();
        return sampleToTimeSpan;
    }
    
    protected Map<String, Integer> loadSampleToSurvivalRate() throws IOException {
        //String fileName = TCGA_GBM_DIR + "TCGASampleSurvivalRateAndClusters.txt";
        String fileName = TCGA_GBM_DIR + "TCGA_Clinical_Info_090909/tcga_GBM_clinical_csv.txt";
        fu.setInput(fileName);
        String line = fu.readLine(); // An empty line
        line = fu.readLine();
        Map<String, Integer> sampleToTime = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            //String[] tokens = line.split("\t");
            String[] tokens = line.split(",");
            if (tokens[0].equals("TCGA-08-0351"))
                continue; // Data is wrong for this sample
            // Last column
            String length = tokens[tokens.length - 1];
            if (length.equals("NOT AVAILABLE"))
                continue;
            sampleToTime.put(tokens[0], new Integer(length));
        }
        fu.close();
        return sampleToTime;
    }
    
    @Test
    public void testLoadSampleToSuvivalRate() throws IOException {
        String fileName = TCGA_GBM_DIR + "TCGASampleSurvivalRateAndClusters.txt";
        // Load samples list
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> samples = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            samples.add(tokens[0]);
        }
        fu.close();
        Map<String, Integer> sampleToSurvival = loadSampleToSurvivalRate();
        Map<String, Integer> sampleToVital = loadSampleToVitalStatus();
        for (String sample : samples) {
            Integer survival = sampleToSurvival.get(sample);
            Integer vital = sampleToVital.get(sample);
            System.out.println(sample + "\t" + survival + "\t" + vital);
        }
    }
    
    /**
     * This method is used to load sample to vital status (living or dead, 
     * 0 for living and 1 for dead).
     * @return
     * @throws IOException
     */
    private Map<String, Integer> loadSampleToVitalStatus() throws IOException {
        String fileName = TCGA_GBM_DIR + "TCGA_Clinical_Info_090909/tcga_GBM_clinical_csv.txt";
        fu.setInput(fileName);
        String line = fu.readLine(); // An empty line
        line = fu.readLine();
        Map<String, Integer> sampleToVitalStatus = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            if (tokens[0].equals("TCGA-08-0351"))
                continue; // Data is wrong for this sample
            // second last column
            String status = tokens[tokens.length - 2].trim();
            if (status.equals("NOT AVAILABLE"))
                continue;
            sampleToVitalStatus.put(tokens[0], new Integer(status));
        }
        fu.close();
        return sampleToVitalStatus;
    }
    
    protected Map<String, String> loadSampleToRecurrent() throws IOException {
        String fileName = TCGA_GBM_DIR + "TCGASampleSurvivalRateAndClusters.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, String> sampleToRecurrent = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            sampleToRecurrent.put(tokens[0], 
                                  tokens[tokens.length - 1]);
        }
        fu.close();
        // Not sure if the following is necessary?
//        // To changes should be made based on the latest version of clinical information from sftp
//        sampleToRecurrent.put("TCGA-02-0014", "No");
//        sampleToRecurrent.put("TCGA-02-0058", "No");
        return sampleToRecurrent;
    }
    
    private List<String> loadHyperMutatedSamples() throws IOException {
        String fileName = TCGA_GBM_DIR + "IndividualSamples.txt";
        fu.setInput(fileName);
        List<String> samples = new ArrayList<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 10)
                continue;
            if (tokens[9].equals("Yes"))
                samples.add(tokens[0]);
        }
        fu.close();
        return samples;
    }
    
    /**
     * This method is used to analyze gene expression data set.
     * @throws Exception
     */
    @Test
    public void clusterGeneDiffGenes() throws Exception {
        Set<String> diffGenes = fu.loadInteractions(TCGA_GBM_DIR + "GBM_EXP_UP_GENES.txt");
        System.out.println("Hyper expression genes: " + diffGenes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        Set<String> fisInGenes = InteractionUtilities.getFIs(diffGenes, fis);
//        fu.saveInteractions(fisInGenes, TCGA_GBM_DIR + "FIsInGBMExpUpGenes.txt");
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String clusterFileName = TCGA_GBM_DIR + "ClustersForExpUpGenes.txt";
        //List<Set<String>> clusters = clusterAnalyzer.cluster(diffGenes);
        //clusterAnalyzer.outputNetworkClusters(clusters, clusterFileName);
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        int index = 0;
        for (Set<String> cluster : clusters) {
            System.out.println(index + "\t" + cluster.size());
            index ++;
        }
        // Check with network clusters
        List<String> samples = loadResequencedSamples();
        Set<String> alteredGenes = getAlteredGenesInSamples(samples);
        System.out.println("Total altered genes: " + alteredGenes.size());
        alteredGenes.retainAll(allGenes);
        System.out.println("Total altered genes in network: " + alteredGenes.size());
        // Calculate FIs
        index = 0;
        System.out.println("\nOverlapping with all altered genes:");
        for (Set<String> cluster : clusters) {
            if (cluster.size() < 5)
                continue;
            Set<String> shared = new HashSet<String>(cluster);
            shared.retainAll(alteredGenes);
            double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(), 
                                                                        cluster.size(),
                                                                        alteredGenes.size(), 
                                                                        shared.size());
            System.out.println(index + "\t" + cluster.size() + "\t" + 
                               shared.size() + "\t" + pvalue + "\t" + 
                               shared);
            index ++;
        }
        System.out.println("\nOverlapping with original network clusters:");
        // Want to get pair-wide network cluster relationships
        List<Set<String>> originalClusters = clusterAnalyzer.loadNetworkClusters(TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster1 = clusters.get(i);
            if (cluster1.size() < 5)
                continue;
            for (int j = 0; j < originalClusters.size(); j++) {
                Set<String> cluster2 = originalClusters.get(j);
                if (cluster2.size() < 5)
                    continue;
                Set<String> shared = new HashSet<String>(cluster1);
                shared.retainAll(cluster2);
                if (shared.size() == 0)
                    continue; // Only want to print out if there is a shared gene
                double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(),
                                                                            cluster1.size(), 
                                                                            cluster2.size(), 
                                                                            shared.size());
                System.out.println(i + "\t" + j + "\t" +
                                   shared.size() + "\t" + pvalue + "\t" + 
                                   shared);
            }
        }
        // Annotate these clusters
        System.out.println("\nPathway annotation for gene exp clusters:");
        clusterAnalyzer.annotateNetworkClusters(clusters, 5);
        System.out.println("\nPathway annotation for original network clusters:");
        clusterAnalyzer.annotateNetworkClusters(originalClusters, 5);
    }
    
    @Test
    public void hierarchicalClusterAlteredGenes() throws Exception {
        // TCGA altered genes
//        List<String> samples = loadResequencedSamples();
//        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
////         Parsons altered genes
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        // Want to have genes have occurring in two or more samples
        Set<String> alteredGenes = selectGenesInSamples(2, sampleToAlteredGenes);
        // To clustering
        // Do a simple sorting to make clustering same since the equal break is actual
        // random
        List<String> geneList = new ArrayList<String>(alteredGenes);
        Collections.sort(geneList);
        hierarchicalClusterGenes(geneList);
    }
    
    private Set<String> loadGBMCore(String dataType) throws IOException {
        String fileName = TCGA_GBM_DIR + "GBMCore.txt";
        Set<String> genes = new HashSet<String>();
        fu.setInput(fileName);
        String line = null;
        boolean inData = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#") && line.contains(dataType)) {
                inData = true;
            }
            else if (inData) {
                String[] tokens = line.split(", ");
                for (String token : tokens)
                    genes.add(token);
                break;
            }
        }
        fu.close();
        return genes;
    }
    
    /**
     * This method is used to check some properties for GBM cores generated from plot
     * based on total genes vs. samples.
     * @throws Exception
     */
    @Test
    public void checkGBMCore() throws Exception {
        // Used to filter out genes not in the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> fisInComp = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genesInFIsInComp = InteractionUtilities.grepIDsFromInteractions(fisInComp);
        
        Set<String> tcgaCore = loadGBMCore("TCGA");
        Set<String> tcgaCoreCopy = new HashSet<String>(tcgaCore);
        System.out.println("Total TCGA Core: " + tcgaCore.size());
        tcgaCore.retainAll(genesInFIs);
        System.out.println("In the FI network: " + tcgaCore.size());
        tcgaCore.retainAll(genesInFIsInComp);
        System.out.println("In the biggest component: " + tcgaCore.size());
        System.out.println("Genes: " + tcgaCore);
        
        Set<String> parsonsCore = loadGBMCore("Parsons");
        Set<String> parsonsCoreCopy = new HashSet<String>(parsonsCore);
        System.out.println("Total Parsons Core: " + parsonsCore.size());
        parsonsCore.retainAll(genesInFIs);
        System.out.println("In the FI network: " + parsonsCore.size());
        parsonsCore.retainAll(genesInFIsInComp);
        System.out.println("In the biggest component: " + parsonsCore.size());
        System.out.println("Genes: " + parsonsCore);
        
        // Check gene sharing
        Set<String> shared = InteractionUtilities.getShared(tcgaCoreCopy, parsonsCoreCopy);
        System.out.println("Total shared genes: " + shared.size());
        System.out.println("Shared genes: " + shared);
        double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES, // Total swiss prot number. This is just an estimate  
                                                                    parsonsCoreCopy.size(), 
                                                                    tcgaCoreCopy.size(),
                                                                    shared.size());
        System.out.println("pvalue using all core genes: " + pvalue);
        // Check on genes in the FI network only
        shared.retainAll(genesInFIsInComp);
        System.out.println("Shared genes in the biggest comp: " + shared.size());
        System.out.println("Shared genes: " + shared);
        pvalue = MathUtilities.calculateHypergeometricPValue(genesInFIsInComp.size(), 
                                                             parsonsCore.size(), 
                                                             tcgaCore.size(), 
                                                             shared.size());
        System.out.println("pvalue using genes in the biggest comp: " + pvalue);
        
//        // Want to calculate span for core
//        String gbmFiFileName = TCGA_GBM_DIR + "FIsForTCGAGBMCore.txt";
//        String gbmSpanFileName = TCGA_GBM_DIR + "SpanForTCGAGBMCore.txt";
//        calculateMinimumSpan(gbmFiFileName, 
//                             gbmSpanFileName,
//                             new ArrayList<String>(tcgaCore));
//        
//        gbmFiFileName = TCGA_GBM_DIR + "FIsForParsonsGBMCore.txt";
//        gbmSpanFileName = TCGA_GBM_DIR + "SpanForParsonsGBMCore.txt";
//        calculateMinimumSpan(gbmFiFileName, 
//                             gbmSpanFileName,
//                             new ArrayList<String>(parsonsCore));
    }
    
    /**
     * This method is used to check total altered genes and average shortest path among altered genes vs 
     * total samples. Genes will be filtered based on the sample numbers.
     * @throws Exception
     */
    @Test
    public void checkAlteredGenesVsSampleCutoff() throws Exception {
//        // TCGA altered genes
//        List<String> samples = loadResequencedSamples();
//        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        // Parsons data set
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes: " + geneToSamples.size());
        Set<String> alteredGenes = new HashSet<String>(geneToSamples.keySet());
        //alteredGenes.retainAll(genesInFIs);
        //System.out.println("In FI network (biggest component): " + alteredGenes.size());
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        List<String> checkingGenes = new ArrayList<String>();
        System.out.println("Samples\tTotal_Genes\tGene_Percentage\tShortestPath\tGenes");
        for (int i = 1; i < sampleToAlteredGenes.size(); i++) {
            // Do a filtering
            checkingGenes.clear();
            for (String gene : alteredGenes) {
                Set<String> samples1 = geneToSamples.get(gene);
                if (samples1.size() >= i)
                    checkingGenes.add(gene);
            }
            if (checkingGenes.size() == 0)
                break;
            //if (checkingGenes.size() < 2)
            //    break;
            // Only working genes in the FI network only
            List<String> checkingGenesInFI = new ArrayList<String>(checkingGenes);
            checkingGenesInFI.retainAll(genesInFIs);
            double shortestPath = 0.0;
            if (checkingGenesInFI.size() > 1)
                shortestPath = calculateShortestPath(checkingGenesInFI,
                                                     bfs,
                                                     nodeToEdges);
            System.out.println(i + "\t" + 
                               checkingGenes.size() + "\t" + 
                               (double)checkingGenes.size() / geneToSamples.size() + "\t" +
                               (shortestPath == 0.0 ? "" : shortestPath) + "\t" +
                               (checkingGenes.size() < 50 ? checkingGenes : ""));
        }
    }
    
    /**
     * This is a permutation test for method checkAlterdGenesVsSampleCutoff.
     * @throws Exception
     */
    @Test
    public void permutationTestAlteredGenesVsSampleCutoff() throws Exception {
        // TCGA altered genes
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
//        // Parsons data set
//        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes: " + geneToSamples.size());
        Set<String> alteredGenes = new HashSet<String>(geneToSamples.keySet());
        //alteredGenes.retainAll(genesInFIs);
        int permutationNumber = 1000;
        // Initialize an integer matrix for altered gene numbers
        List<int[]> sampleToGeneNumbers = new ArrayList<int[]>();
        for (int i = 0; i < sampleToAlteredGenes.size(); i++)
            sampleToGeneNumbers.add(new int[permutationNumber]);
        // Initialize an double matrix for average shortest path
        List<double[]> sampleToAvgPaths = new ArrayList<double[]>();
        for (int i = 0; i < sampleToAlteredGenes.size(); i++)
            sampleToAvgPaths.add(new double[permutationNumber]);
        for (int perm = 0; perm < permutationNumber; perm++) {
            System.out.println("Permutation test: " + perm);
            // Random sampling
            //Set<String> randomGenes = MathUtilities.randomSampling(genesInFIs, alteredGenes.size());
            // Re-sampling
            Map<String, Set<String>> randomSampleToAlteredGenes = new HashMap<String, Set<String>>();
            for (String sample : sampleToAlteredGenes.keySet()) {
                Set<String> set = sampleToAlteredGenes.get(sample);
                Set<String> sampled = MathUtilities.randomSampling(alteredGenes, 
                                                                   set.size());
                randomSampleToAlteredGenes.put(sample, sampled);
            }
            Map<String, Set<String>> randomGeneToSamples = InteractionUtilities.switchKeyValues(randomSampleToAlteredGenes);
            List<String> checkingGenes = new ArrayList<String>();          
            for (int i = 1; i < randomSampleToAlteredGenes.size(); i++) {
                // Do a filtering
                checkingGenes.clear();
                for (String gene : randomGeneToSamples.keySet()) {
                    Set<String> samples1 = randomGeneToSamples.get(gene);
                    if (samples1.size() >= i)
                        checkingGenes.add(gene);
                }
                if (checkingGenes.size() == 0)
                    break;
                double shortestPath = 0.0d;
                List<String> checkingGenesInFIs = new ArrayList<String>(checkingGenes);
                checkingGenesInFIs.retainAll(genesInFIs);
                if (checkingGenesInFIs.size() > 1)
                    shortestPath = calculateShortestPath(checkingGenesInFIs,
                                                         bfs,
                                                         nodeToEdges);
                int[] geneNumbers = sampleToGeneNumbers.get(i - 1);
                geneNumbers[perm] = checkingGenes.size();
                double[] avgPaths = sampleToAvgPaths.get(i - 1);
                avgPaths[perm] = shortestPath;
            }
        }
        String outFileName = TCGA_GBM_DIR + "PermTestAlteredGenesVsSamples_TCGA_1000_No_0.txt";
        fu.setOutput(outFileName);
        // Generate output
        DescriptiveStatistics stat = new DescriptiveStatistics();
        fu.printLine("Sample_Cutoff\tGene_Number\tSD_Gene_Number\tShortest_Path\tSD_Shortest_Path");
        int[] totalGenes = sampleToGeneNumbers.get(0);
        for (int i = 0; i < sampleToGeneNumbers.size(); i++) {
            int[] geneNumbers = sampleToGeneNumbers.get(i);
            // Want to convert it to percentage
            for (int j = 0; j < geneNumbers.length; j++) {
                stat.addValue((double) geneNumbers[j] / totalGenes[j]);
            }
            double avgNumberPercentage = stat.getMean();
            double sdAvgNoPer = stat.getStandardDeviation();
            double[] avgPaths = sampleToAvgPaths.get(i);
            stat.clear();
            for (int j = 0; j < avgPaths.length; j++) {
                // Don't count zero for average path since it is meaingless
                if (avgPaths[j] > 0.0d)
                    stat.addValue(avgPaths[j]);
            }
            double avgPath = stat.getMean();
            double sdAvgPath = stat.getStandardDeviation();
            fu.printLine((i + 1) + "\t" + 
                         avgNumberPercentage + "\t" +
                         sdAvgNoPer + "\t" +
                         avgPath + "\t" +
                         sdAvgPath);
            stat.clear();
        }
        fu.close();
    }
    
    @Test
    public void clusterMutatedGenes() throws Exception {
//        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        NetworkClusterAnalyzer networkAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> mutatedClusters = networkAnalyzer.loadNetworkClusters(TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes.txt");
//        List<Set<String>> mutatedClusters = networkAnalyzer.cluster(mutatedGenes, 0.2);
//        for (int i = 0; i < mutatedClusters.size(); i++) {
//            Set<String> cluster = mutatedClusters.get(i);
//            System.out.println(i + "\t" + cluster.size());
//        }
//        String outFileName = TCGA_GBM_DIR + "MutatedGenes2OrMoreClusters.txt";
//        networkAnalyzer.outputNetworkClusters(mutatedClusters, 
//                                              outFileName);
        // Generate a a list of FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        Set<String> diffFIs = InteractionUtilities.getFIs(mutatedGenes, fis);
//        outFileName = TCGA_GBM_DIR + "FIsAmongMutatedGenes2OrMore.txt";
//        fu.saveInteractions(diffFIs, outFileName);
        // Check distance between two clusters
        List<Set<String>> expClusters = networkAnalyzer.loadNetworkClusters(TCGA_GBM_DIR + "ExpUpDownGenesNetworkClusters.txt");
        // Check distance between two clusters
        System.out.println("\nDistance from mutated clusters to gene-up clusters:");
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(fis);
        for (int i = 0; i < mutatedClusters.size(); i++) {
            Set<String> mutatedCluster = mutatedClusters.get(i);
            if (mutatedCluster.size() < 5)
                continue;
            for (int j = 0; j < expClusters.size(); j++) {
                Set<String> geneUpCluster = expClusters.get(j);
                if (geneUpCluster.size() < 5)
                    continue;
                double dist = calculateAverageDistance(mutatedCluster, 
                                                       geneUpCluster, 
                                                       bfs,
                                                       idToPartners);
                System.out.println(i + "-" + j + "\t" + dist);
            }
        }
    }
    
    /**
     * Analyze properties differentially expressed genes based on the FIs.
     */
    @Test
    public void analyzeDiffExpGenes() throws Exception {
        Collection<String> diffGenes = loadExpDiffGenes();
        NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        List<Set<String>> diffClusters = analyzer.cluster(diffGenes, 0.1);
        for (int i = 0; i < diffClusters.size(); i++) {
            System.out.println(i + "\t" + diffClusters.get(i).size());
        }
        // Output
        //String outFileName = TCGA_GBM_DIR + "ExpUpGenesNetworkClusters.txt";
        //String outFileName = TCGA_GBM_DIR + "ExpDownGenesNetworkClusters.txt";
        String outFileName = TCGA_GBM_DIR + "ExpUpDownGenesNetworkClusters.txt";
        analyzer.outputNetworkClusters(diffClusters,
                                       outFileName);
        // Generate a a list of FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        diffGenes.retainAll(fiGenes);
        System.out.println("diff genes in FIs: " + diffGenes.size());
        Set<String> diffFIs = InteractionUtilities.getFIs(diffGenes, fis);
        // outFileName = TCGA_GBM_DIR + "FIsInExpUpGenes.txt";
        // outFileName = TCGA_GBM_DIR + "FIsInExpDownGenes.txt";
        outFileName = TCGA_GBM_DIR + "FIsInExpUpDownGenes.txt";
        fu.saveInteractions(diffFIs, outFileName);
        // Check graph components
        List<Set<String>> components = new GraphAnalyzer().calculateGraphComponents(diffFIs);
        // print out
        System.out.println("\nConnected graph component: " + components.size());
        for (int i = 0; i < components.size(); i++) {
            System.out.println(i + "\t" + components.get(i).size());
        }
        // Check average shortest path for the biggest component
        Set<String> biggest = components.get(0);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        double avgPath = calculateShortestPath(new ArrayList<String>(biggest),
                                               bfs,
                                               nodeToEdges);
        System.out.println("Avg shortest path: " + avgPath);
    }
    
    /**
     * This method is used to do a hierarchically clustering on the diff genes.
     * @throws IOException
     */
    @Test
    public void hierarchicalClusterExpDiffGenes() throws IOException {
        Collection<String> diffGenes = loadExpDiffGenes();
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        List<String> geneList = new ArrayList<String>(diffGenes);
        System.out.println("Total diff genes: " + geneList.size());
        geneList.retainAll(geneToPartners.keySet());
        System.out.println("Diff genes in the FI network: " + geneList.size());
        Map<String, Integer> pairToDistance = calculateDistances(geneList, 
                                                                 bfs, 
                                                                 geneToPartners);
        List<HierarchicalClusterNode> clusters = hierarchicalCluster(geneList, 
                                                         pairToDistance,
                                                         false);
        HierarchicalClusterNode root = clusters.get(0);
        outputCluster(root);
    }
    
    protected Collection<String> loadExpDiffGenes() throws IOException {
        String fileName = TCGA_GBM_DIR + "Gene_exp_p_t_values.txt";
        final Map<String, Double> geneToTScore = new HashMap<String, Double>();
        Map<String, Double> genetoPValue = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double tscore = new Double(tokens[1]);
            geneToTScore.put(tokens[0], tscore);
            genetoPValue.put(tokens[0], new Double(tokens[2]));
        }
        fu.close();
        // Sort based on t-score
        List<String> geneList = new ArrayList<String>(geneToTScore.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double tscore1 = geneToTScore.get(gene1);
                Double tscore2 = geneToTScore.get(gene2);
                return tscore2.compareTo(tscore1);
            }
        });
//        int c = 1;
//        for (String gene : geneList) {
//            System.out.println(c + "\t" + gene + "\t" + geneToTScore.get(gene) + "\t" + genetoPValue.get(gene));
//            c ++;
//        }
        // Pick up the first two hundred genes
        //int total = 400;
        // Use cutoff value p-value = 1.0e-6
        Double pvalueCutoff = 1.0E-6;
        List<String> diffGenes = new ArrayList<String>();
        for (String gene : geneList) {
            Double tscore = geneToTScore.get(gene);
            if (tscore > 0) {// For up genes
                Double pvalue = genetoPValue.get(gene);
                if (pvalue < pvalueCutoff)
                    diffGenes.add(gene);
            }
            else if (tscore < 0) { // For down genes
                Double pvalue = genetoPValue.get(gene);
                if (pvalue < pvalueCutoff)
                    diffGenes.add(gene);
            }
        }
        System.out.println("Total diff genes: " + diffGenes.size());
        return diffGenes;
    }
    
    /**
     * This method is used to do permutation test to prove that 44% of diff genes
     * can form a single connected graph component is not by random.
     * @throws Exception
     */
    @Test
    public void permutateGraphComponentForArrayGenes() throws Exception {
        // To load array genes
        String fileName = TCGA_GBM_DIR + "Gene_exp_p_t_values.txt";
        fu.setInput(fileName);
        Set<String> arrayGenes = new HashSet<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            arrayGenes.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total array genes: " + arrayGenes.size());
        // FI genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        arrayGenes.retainAll(fiGenes);
        System.out.println("Array genes in FI network: " + arrayGenes.size());
        // Starting permutation
        int permutation = 1000;
        int setSize = 462; // Diff genes in the FI network
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Double> values = new ArrayList<Double>();
        //BreadthFirstSearch bfs = new BreadthFirstSearch();
        //Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        for (int i = 0; i < permutation; i++) {
            Set<String> sample = MathUtilities.randomSampling(arrayGenes, setSize);
            Set<String> interactions = InteractionUtilities.getFIs(sample, fis);
            List<Set<String>> components = graphAnalyzer.calculateGraphComponents(interactions);
            double value = (double) components.get(0).size() / sample.size();
            values.add(value);
//            double path = calculateShortestPath(new ArrayList<String>(components.get(0)),
//                                                bfs,
//                                                nodeToEdges);
//            System.out.println(value + "\t" + path);
        }
        Collections.sort(values, new Comparator<Double>() {
            public int compare(Double value1, Double value2) {
                return value2.compareTo(value1);
            }
        });
        for (int i = 0; i < values.size(); i++) {
            System.out.println((i + 1) + "\t" + values.get(i));
        }
    }
    
    /**
     * This method is used to calculate centralities of genes in a network cluster.
     * @throws Exception
     */
    @Test
    public void calculateCentralityInNetworkClusters() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        //List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes.txt");
        List<Set<String>> clusters = new MCLResultsAnalyzer().loadMCLClusters(R3Constants.RESULT_DIR + "MCLCluster_FIsInGene_041709_I40.txt");
        NetworkClusterAnalyzer networkAnalyzer = new NetworkClusterAnalyzer();
        JungGraphUtilities gu = new JungGraphUtilities();
        int index = 0;
        for (Set<String> cluster : clusters) {
            System.out.println("Cluster " + index);
            index ++;
            Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
            final Map<String, Double> geneToCentrality = gu.calculateCentralities(cluster, fisInCluster);
            List<String> geneList = new ArrayList<String>(geneToCentrality.keySet());
            Collections.sort(geneList, new Comparator<String>() {
                public int compare(String gene1, String gene2) {
                    Double value1 = geneToCentrality.get(gene1);
                    Double value2 = geneToCentrality.get(gene2);
                    return value2.compareTo(value1);
                }
            });
            for (String gene : geneList) {
                Double centrality = geneToCentrality.get(gene);
                if (centrality >= 0.005)
                    System.out.println(gene + "\t" + geneToCentrality.get(gene));
            }
            System.out.println();
        }
    }
    
    @Test
    public void annotateNetworkClusters() throws Exception {
        //String clusterFile = TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes.txt";
        //String clusterFile = TCGA_GBM_DIR + "ExpUpDownGenesNetworkClusters.txt";
        //String  clusterFile = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        //String clusterFile = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        String clusterFile = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesInPathwayFIs080509_220.txt";
        NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = analyzer.loadNetworkClusters(clusterFile);
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            System.out.println(i + "\t" + cluster.size());
        }
        analyzer.annotateNetworkClusters(clusterFile, 5);
//        Map<String, List<Set<String>>> cutoffToNetworkClusters = loadMiniNetworkClusters();
//        for (String cutoff : cutoffToNetworkClusters.keySet()) {
//            if (cutoff.contains("6")) {
//                List<Set<String>> clusters = cutoffToNetworkClusters.get(cutoff);
//                analyzer.annotateNetworkClusters(clusters, 5);
//            }
//        }
    }
    
    @Test
    public void annotataeNetworkClustersWithGO() throws Exception {
        //String clusterFile = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        //NetworkClusterAnalyzer analyzer = new NetworkClusterAnalyzer();
        //List<Set<String>> clusters = analyzer.loadNetworkClusters(clusterFile);
        List<Set<String>> clusters = null;
        Map<String, List<Set<String>>> cutoffToNetworkClusters = loadMiniNetworkClusters();
        for (String cutoff : cutoffToNetworkClusters.keySet()) {
            if (cutoff.contains("6")) {
                clusters = cutoffToNetworkClusters.get(cutoff);
            }
        }
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        int index = -1;
        for (Set<String> cluster : clusters) {
            index ++;
            if (cluster.size() < 5)
                continue;
            System.out.println("\nCluster " + index + ": " + cluster.size());
            annotator.annotateGenesWithFDR(cluster, AnnotationType.BP);
        }
    }
    
    private List<Set<String>> generateRandomClusters(List<Set<String>> clusters,
                                                     Set<String> totalGenes,
                                                     RandomData randomizer) {
        List<Set<String>> samples = new ArrayList<Set<String>>();
        Set<String> copy = new HashSet<String>(totalGenes);
        for (Set<String> cluster : clusters) {
            int size = cluster.size();
            Set<String> tmp = MathUtilities.randomSampling(copy, 
                                                           size,
                                                           randomizer);
            samples.add(tmp);
            copy.removeAll(tmp);
        }
        return samples;
    }
    
    /**
     * In this permutation test, the clusters genes are randomized.
     * @throws Exception
     */
    @Test
    public void permutatitonTestForSampleInNetworkClusters2() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        // To test clusters from Parsons', random genes should use genes from parson's data set.
        sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        // Test for Science GBM genes
//        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();        
        // Will focus on FI genes only since only we know these genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> totalGenes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            genes.retainAll(fiGenes);
        }
        for (Set<String> set : sampleToAlteredGenes.values())
            totalGenes.addAll(set);
        System.out.println("Total genes: " + totalGenes.size());
        String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(clusterFileName);
        // Test for Science GBM genes
        // For science permutation, the total genes should use TCGA altered genes in FI since
        // the null hypothesis is that the first cluster0 and 1 are not significant.
        //sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        // The target data set for test
        sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        double target = 0.7142857142857143;
        int permutationNumber = 10000;
        List<Double> permutPercentages = new ArrayList<Double>();
        RandomData randomizer = new RandomDataImpl();
        for (int j = 0; j < permutationNumber; j++) {
            List<Set<String>> randomClusters = generateRandomClusters(clusters, 
                                                                      totalGenes, 
                                                                      randomizer);
            int counter = 0;
            for (String sample : sampleToAlteredGenes.keySet()) {
                Set<String> genes = sampleToAlteredGenes.get(sample);
                // Want to get the list of clusters
                List<Integer> clusterIds = new ArrayList<Integer>();
                for (int i = 0; i < randomClusters.size(); i++) {
                    Set<String> set = randomClusters.get(i);
                    for (String gene : genes) {
                        if (set.contains(gene)) {
                            clusterIds.add(i);
                            break;
                        }
                    }
                }
                if (isSampleInClusters(clusterIds))
                    counter ++;
            }
            double percent = (double) counter / sampleToAlteredGenes.size();
            permutPercentages.add(percent);
        }
        Collections.sort(permutPercentages, new Comparator<Double>() {
            public int compare(Double v1, Double v2) {
                return v2.compareTo(v1);
            }
        });
        //for (int i = 0; i < permutPercentages.size(); i++)
        //    System.out.println(i + ": " + permutPercentages.get(i));
        calculatePValue(target, 
                        permutationNumber,
                        permutPercentages);
    }
    
    @Test
    public void permutationTestForSamplesInNetworkClusters() throws Exception {
        // To test resequenced samples for TCGA
//        List<String> samples = loadResequencedSamples();
//        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
//        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
//        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToCNVGenes(samples);
        // For testing Science GBM genes
        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        summarizeSampleToGenes(sampleToAlteredGenes);
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesInPathwayFIs080509_115.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes042810.txt";
        double target = 0.7142857142857143;
//        double target = 0.8351648351648352;
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
//        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
//        for (Set<String> cluster : clusters) {
//            System.out.println(cluster.size());
//        }
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);;
        clusterAnalyzer.permutationTestForSamplesInNetworkClusters(sampleToAlteredGenes,
                                                                   clusterFileName, 
                                                                   targetClusters, 
                                                                   target,
                                                                   false);
    }
    
    /**
     * Permutation tests considering edge betweenness.
     * @throws Exception
     */
    @Test
    public void permutationTestForSamplesInNetworkClustersOnEdgeBetweenness() throws Exception {
        // To test resequenced samples for TCGA
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        //        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        //        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToCNVGenes(samples);
        // For testing Science GBM genes
//        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
//        summarizeSampleToGenes(sampleToAlteredGenes);
        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesInPathwayFIs080509_115.txt";
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        //        double target = 0.7142857142857143;
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        for (Set<String> cluster : clusters) {
            System.out.println(cluster.size());
        }
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        fisGenes.retainAll(InteractionUtilities.grepAllGenes(sampleToAlteredGenes));
        System.out.println("Total genes: " + fisGenes.size());
        //        Set<String> fisInGenes = InteractionUtilities.getFIs(InteractionUtilities.grepAllGenes(sampleToAlteredGenes),
        //                                                             fis);
        //        double modularity = clusterAnalyzer.calculateModualarity(clusters, fisInGenes);
        //        System.out.println("Modularity: " + modularity);
        //        clusterAnalyzer.permutationTestBasedOnModularity(sampleToAlteredGenes);
        List<Integer> targetClusters = new ArrayList<Integer>();
        targetClusters.add(0);
        targetClusters.add(1);;
        List<Integer> targetClusterSizes = new ArrayList<Integer>();
        for (Integer index : targetClusters) {
            Set<String> cluster = clusters.get(index);
            targetClusterSizes.add(cluster.size());
        }
        clusterAnalyzer.permutationTestForSamplesInClustersOnBetweenness(sampleToAlteredGenes, 
                                                                         targetClusters,
                                                                         targetClusterSizes);
    }

    private void calculatePValue(double target, 
                                 int permutationNumber,
                                 List<Double> permutPercentages) {
        int index = -1;
        for (int i = 0; i < permutPercentages.size(); i++) {
            double pert = permutPercentages.get(i);
            if (pert < target) {
                index = i;
                break;
            }
        }
        if (index == 0)
            System.out.println("pvalue < " + 1.0 / permutationNumber);
        else
            System.out.println("pvalue: " + (double) index / permutationNumber);
        double max = permutPercentages.get(0);
        System.out.println("Max value: " + max);
        System.out.println("Target value: " + target);
    }
    
    private boolean isSampleInClusters(List<Integer> clusterIds) {
        // Check for a sample has mutation in both cluster 0 and 1.
        if (clusterIds.contains(0) && 
            clusterIds.contains(2))// && 
//            clusterIds.contains(2))
            return true;
//         Check for a sample has a mutation in either cluster 0 or 1.
//        if (clusterIds.contains(0) || clusterIds.contains(1))
//            return true;
//        if (clusterIds.contains(0))
//            return true;
        return false;
    }
    
    /**
     * This method is used to check altered genes in some samples.
     * @throws Exception
     */
    @Test
    public void checkSampleAlteredGenes() throws Exception {
        String sampleNames = "TCGA-06-0201, TCGA-02-0047, TCGA-06-0219, TCGA-02-0024, " +
        		"TCGA-02-0038, TCGA-06-0143, TCGA-06-0154, TCGA-02-0054, " +
        		"TCGA-02-0115, TCGA-06-0122, TCGA-06-0158, TCGA-06-0166, " +
        		"TCGA-06-0133, TCGA-02-0027, TCGA-02-0086, TCGA-06-0241, " +
        		"TCGA-06-0210, TCGA-02-0116, TCGA-06-0171, TCGA-06-0185, " +
        		"TCGA-06-0156, TCGA-02-0034, TCGA-06-0126";
        String[] samples = sampleNames.split(", ");
        System.out.println("Total samples: " + samples.length);
        Map<String, Integer> sampleToSurvival = loadSampleToTimeSpan();
        List<String> clusterSamples = new ArrayList<String>();
        for (String sample : samples) {
            Integer survival = sampleToSurvival.get(sample);
            //System.out.println(sample + "\t" + survival);
            if (survival < 2)
                clusterSamples.add(sample);
        }
        List<String> resequencedSamples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(resequencedSamples);
        Set<String> clusterSampleGenes = new HashSet<String>();
        for (String sample : clusterSamples) {
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            clusterSampleGenes.addAll(alteredGenes);
        }
        System.out.println("Altered genes in short samples: " + clusterSampleGenes);
        // Get genes for long short samples
        List<String> longSampleList = new ArrayList<String>();
        List<String> shortSampleList = new ArrayList<String>();
        for (String sample : resequencedSamples) {
            Integer survival = sampleToSurvival.get(sample);
            if (survival < 2)
                shortSampleList.add(sample);
            else if (survival >= 2)
                longSampleList.add(sample);
        }
        System.out.println("Short samples: " + shortSampleList.size());
        System.out.println("Long samples: " + longSampleList.size());
        Set<String> longSampleGenes = new HashSet<String>();
        for (String sample : longSampleList) {
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            longSampleGenes.addAll(alteredGenes);
        }
        System.out.println("Long sample genes: " + longSampleGenes);
        // Calculate distance to this core from the cluster sample genes
        System.out.println("Short samples:");
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        clusterSampleGenes.retainAll(fiGenes);
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(fis);
        for (String sample : shortSampleList) {
            if (clusterSamples.contains(sample))
                continue;
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            alteredGenes.retainAll(fiGenes);
            double dist = calculateMinShortestPathOneWay(alteredGenes,
                                                         clusterSampleGenes,
                                                         bfs, 
                                                         idToPartners);
            System.out.println(sample + ": " + dist);
        }
        System.out.println("Long samples:");
        for (String sample : longSampleList) {
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            if (alteredGenes == null)
                continue;
            alteredGenes.retainAll(fiGenes);
            double dist = calculateMinShortestPathOneWay(alteredGenes,
                                                         clusterSampleGenes,
                                                         bfs, 
                                                         idToPartners);
            System.out.println(sample + ": " + dist);
        }
    }
    
    /**
     * This method is used to analyze altered genes from samples and try to find a minimal set of
     * genes for both cluster 0 and 1 in order to figure out the minimal requirement of alteration
     * for a cancer sample.
     * @throws Exception
     */
    @Test
    public void checkSampleAlteredGenesInNetworkClusters() throws Exception {
        // TCGA data set
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
//        // Parsons data set
//        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
        List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(clusterFileName);

        
        System.out.println("Total clusters: " + clusters.size());
        // Get the list of genes for cluster 0 and 1
        final Map<String, Integer> cluster0GeneToSamples = new HashMap<String, Integer>();
        final Map<String, Integer> cluster1GeneToSamples = new HashMap<String, Integer>();
        Set<String> cluster0 = clusters.get(0);
        Set<String> cluster1 = clusters.get(2);
        // Check how many samples having 0 and 1 clusters
        int counter = 0;
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            for (String gene : genes) {
                if (cluster0.contains(gene)) {
                    Integer c = cluster0GeneToSamples.get(gene);
                    if (c == null)
                        cluster0GeneToSamples.put(gene, 1);
                    else
                        cluster0GeneToSamples.put(gene, ++c);
                }
                if (cluster1.contains(gene)) {
                    Integer c = cluster1GeneToSamples.get(gene);
                    if (c == null)
                        cluster1GeneToSamples.put(gene, 1);
                    else
                        cluster1GeneToSamples.put(gene, ++c);
                }
            }
        }
        // Check total genes in cluster 0 and 1 touched by samples
        System.out.println("Sample touched genes in cluster 0: " + cluster0GeneToSamples.size());
        List<String> geneList0 = outputGeneToSampleCounts(cluster0GeneToSamples);
        System.out.println("Sample touched genes in cluster 1: " + cluster1GeneToSamples.size());
        List<String> geneList1 = outputGeneToSampleCounts(cluster1GeneToSamples);
        // The target
        double targetPercent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                                cluster0,
                                                                cluster1);
        System.out.println("Target percentage: " + targetPercent);
        Set<String> newCluster0 = new HashSet<String>(cluster0);
        Set<String> newCluster1 = new HashSet<String>(cluster1);
        for (int i = geneList0.size() - 1; i >= 0; i--) {
            String gene = geneList0.get(i);
            newCluster0.remove(gene);
            double newPercent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                                 newCluster0,
                                                                 newCluster1);
            if (newPercent < targetPercent)
                break;
        }
        for (int i = geneList1.size() - 1; i >= 0; i--) {
            String gene = geneList1.get(i);
            newCluster1.remove(gene);
            double newPercent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                                 newCluster0,
                                                                 newCluster1);
            if (newPercent < targetPercent)
                break;
        }
        System.out.println("New clusters: ");
        System.out.println("Cluster 0: " + newCluster0.size() + " " + newCluster0);
        System.out.println("Cluster 1: " + newCluster1.size() + " " + newCluster1);
        // Create two new clusters with samples three or more
        //int sampleCutOff = 4;
        for (int sampleCutOff = 1; sampleCutOff < 10; sampleCutOff ++) {
            System.out.println("Cutoff: " + sampleCutOff);
            newCluster0.clear();
            newCluster1.clear();
            for (String gene : cluster0GeneToSamples.keySet()) {
                Integer c = cluster0GeneToSamples.get(gene);
                if (c >= sampleCutOff)
                    newCluster0.add(gene);
            }
            for (String gene : cluster1GeneToSamples.keySet()) {
                Integer c = cluster1GeneToSamples.get(gene);
                if (c >= sampleCutOff)
                    newCluster1.add(gene);
            }
            double newPercent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                                 newCluster0,
                                                                 newCluster1);
            System.out.println("New clusters: ");
            System.out.println("Cluster 0: " + newCluster0.size() + " " + newCluster0);
            System.out.println("Cluster 1: " + newCluster1.size() + " " + newCluster1);
            System.out.println("New percentage: " + newPercent + "\n");
        }
    }
    
    private double calculatePercentageInTwoClusters(Map<String, Set<String>> sampleToAlteredGenes,
                                                    Set<String> cluster0,
                                                    Set<String> cluster1) {
        int hitSamples = 0;
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> altered = sampleToAlteredGenes.get(sample);
            Set<String> shared0 = new HashSet<String>(altered);
            shared0.retainAll(cluster0);
            Set<String> shared1 = new HashSet<String>(altered);
            shared1.retainAll(cluster1);
            //System.out.println(sample + "\t" + shared0 + "\t" + shared1);
            if (shared0.size() > 0 && shared1.size() > 0)
                hitSamples ++;
        }
        return (double) hitSamples / sampleToAlteredGenes.size();
    }
    
    private List<String> outputGeneToSampleCounts(final Map<String, Integer> geneToSamples) {
        List<String> geneList = new ArrayList<String>(geneToSamples.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Integer c1 = geneToSamples.get(gene1);
                Integer c2 = geneToSamples.get(gene2);
                return c2.compareTo(c1);
            }
        });
        for (String gene : geneList)
            System.out.println(gene + "\t" + geneToSamples.get(gene));
        return geneList;
    }
    
    private Map<String, List<Set<String>>> loadMiniNetworkClusters() throws IOException {
        String clusterFileName = TCGA_GBM_DIR + "MiniCluster0And1.txt";
        fu.setInput(clusterFileName);
        String line = null;
        Map<String, List<Set<String>>> cutoffToClusters = new HashMap<String, List<Set<String>>>();
        String cutoffLine = null;
        List<Set<String>> clusters = new ArrayList<Set<String>>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Cutoff")) {
                cutoffLine = line;
                clusters = new ArrayList<Set<String>>();
                cutoffToClusters.put(cutoffLine, clusters);
            }
            else if (line.startsWith("Cluster")) {
                int index = line.indexOf("[");
                String genes = line.substring(index, line.length() - 1);
                String[] tokens = genes.split(", ");
                Set<String> cluster = new HashSet<String>();
                for (String token : tokens)
                    cluster.add(token);
                clusters.add(cluster);
            }
        }
        fu.close();
        return cutoffToClusters;
    }
    
    /**
     * This method is used to check sample distributions for clusters with sample cutoff
     * generated.
     * @throws Exception
     */
    @Test
    public void checkSamplesInMiniNetworkClusters() throws Exception {
        Map<String, List<Set<String>>> cutoffToClusters = loadMiniNetworkClusters();
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        System.out.println("Check TCGA GBM samples: " + sampleToAlteredGenes.size());
        List<Set<String>> clusters;
        for (String cutoff : cutoffToClusters.keySet()) {
            if (!cutoff.contains("6"))
                continue;
            clusters = cutoffToClusters.get(cutoff);
            double percent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                              clusters.get(0),
                                                              clusters.get(1));
            System.out.println(cutoff + "\t" + percent);
        }
        // Check Parsons data sample
        // Altered genes for 22 discovery genes only.
        sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();
        System.out.println("\nCheck Parsons GBM samples: " + sampleToAlteredGenes.size());
        for (String cutoff : cutoffToClusters.keySet()) {
            if (!cutoff.contains("6"))
                continue;
            clusters = cutoffToClusters.get(cutoff);
            double percent = calculatePercentageInTwoClusters(sampleToAlteredGenes, 
                                                              clusters.get(0),
                                                              clusters.get(1));
            System.out.println(cutoff + "\t" + percent);
        }
    }
    
    /**
     * Check how many samples are covered for each clusters from hierarchical clustering.
     * @throws Exception
     */
    @Test
    public void checkSamplesInHierarchicalClusters() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        String clusterFileName = TCGA_GBM_DIR + "HierarchicalClusterResultsForTCGA092909.txt";
        HierarchicalCluster clusterizer = new HierarchicalCluster();
        List<HierarchicalClusterNode> clusterNodes = clusterizer.loadHierarchicalClusters(clusterFileName);
        // Search the samples in clusters
        for (HierarchicalClusterNode node : clusterNodes) {
            // Get genes
            int touched = 0;
            for (String sample : sampleToAlteredGenes.keySet()) {
                Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
                Set<String> shared = InteractionUtilities.getShared(alteredGenes, node.ids);
                if (shared.size() > 1)
                    touched ++;
            }
            System.out.println(node.pathDistance + ": " + touched);
        }
    }
    
    /**
     * This method is used to draw a heat map for a list of samples ordered based on
     * a hierarchical clustering.
     * @throws Exeption
     */
    @Test
    public void drawSampleListInHierarchicalClustersFromR() throws Exception {
        //String fileName = TCGA_GBM_DIR + "SamplesOrderByHclusterFromR.txt";
        String fileName = TCGA_GBM_DIR + "SamplesOrderByHclusterFromRFromAvg.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        fu.close();
        String[] samples = line.split(" ");
        System.out.println("Total samples: " + samples.length);
        Map<String, String> sampleToRec = loadSampleToRecurrent();
        int length = 1116;
        int width  = 50;
        List<String> sampleList = new ArrayList<String>();
        for (String sample : samples)
            sampleList.add(sample);
        HierarchicalCluster drawer = new HierarchicalCluster();
        List<String> infoList = new ArrayList<String>();
        infoList.add("No");
        infoList.add("Rec");
        infoList.add("Sec");
        List<Color> colorList = new ArrayList<Color>();
        colorList.add(Color.green);
        colorList.add(Color.blue);
        colorList.add(Color.red);
        String outFileName = TCGA_GBM_DIR + "SampleListShadedOnRecFromAvg.png";
        drawer.drawSampleInformation(sampleToRec,
                                     sampleList,
                                     infoList,
                                     colorList, 
                                     length,
                                     width, 
                                     outFileName);
    }
    
    /**
     * This method is used to check some features for sample clusters generated from R functio hCluster
     * based on heatmap. The clustering is based on binary distance from network clusters.
     * @throws Exception
     */
    @Test
    public void checkSamplesInHierarchicalClustersFromR() throws Exception {
       String fileName = TCGA_GBM_DIR + "SamplesOrderByHclusterFromR.txt";
       fu.setInput(fileName);
       String line = fu.readLine();
       fu.close();
       String[] samples = line.split(" ");
       System.out.println("Total samples: " + samples.length);
       Map<String, String> sampleToRec = loadSampleToRecurrent();
       int totalNoRec = 0;
       for (String sample : samples) {
           String rec = sampleToRec.get(sample);
           System.out.println(sample + "\t" + rec);
           if (rec.equals("No"))
               totalNoRec ++;
       }
       // Check recurrence label
       String[] bounds = new String[] {
               "TCGA-02-0024",// "TCGA-02-0055",
               "TCGA-06-0154",// "TCGA-02-0114",
               "TCGA-06-0241",// "TCGA-06-0141",
               "TCGA-06-0147",// "TCGA-06-0124",
               "TCGA-06-0139",// "TCGA-06-0189"
       };
       List<String> boundsList = Arrays.asList(bounds);
       List<List<String>> sampleClusters = new ArrayList<List<String>>();
       // Fill cluster
       for (int i = 0; i < bounds.length; i++)
           sampleClusters.add(new ArrayList<String>());
       List<String> cluster = null;
       for (String sample : samples) {
           int index = boundsList.indexOf(sample);
           if (index >= 0) {
               cluster = sampleClusters.get(index);
           }
           cluster.add(sample);
       }
       // Check if the clusters is correct
       for (List<String> sampleCluster : sampleClusters) {
           System.out.println(sampleCluster.size() + ": " + sampleCluster);
       }
       // Do a fisher test
       FisherExact fisher = new FisherExact(100);
       int index = 0;
       for (List<String> sampleCluster : sampleClusters) {
           System.out.printf("%nCluster %d: %d%n",
                             index ++,
                             sampleCluster.size());
           List<String> noInCluster = new ArrayList<String>();
           for (String sample : sampleCluster) {
               String rec = sampleToRec.get(sample);
               if (rec.equals("No"))
                   noInCluster.add(sample);
           }
           int a = noInCluster.size();
           int b = sampleCluster.size() - a;
           int c = totalNoRec - a;
           int d = samples.length - sampleCluster.size() - c;
           double pvalue = fisher.getTwoTailedP(a,
                                                b,
                                                c,
                                                d);
           System.out.printf("Contingence table: %d, %d, %d, %d%n", a, b, c, d);
           System.out.println("pvalue from Fisher: " + pvalue);
           pvalue = MathUtilities.calculateHypergeometricPValue(samples.length, 
                                                                sampleCluster.size(), 
                                                                totalNoRec, 
                                                                noInCluster.size());
           if (pvalue > 0.50) // Try to see another way
               pvalue = MathUtilities.calculateHypergeometricPValue(samples.length, 
                                                                    sampleCluster.size(), 
                                                                    samples.length - totalNoRec, 
                                                                    sampleCluster.size() - noInCluster.size());
           System.out.println("pvalue from hyper: " + pvalue);
       }
       // Check the distribution of hypermutated samples
       List<String> hyperMutatedSamples = loadHyperMutatedSamples();
       System.out.println("\nTotal hypermutated samples: " + hyperMutatedSamples.size());
       index = 0;
       for (List<String> sampleCluster : sampleClusters) {
           System.out.printf("%nCluster %d: %d%n",
                             index ++,
                             sampleCluster.size());
           List<String> hyperInCluster = new ArrayList<String>(sampleCluster);
           hyperInCluster.retainAll(hyperMutatedSamples);
           System.out.println("Hyper-mutated samples in cluster: " + hyperInCluster);
           int a = hyperInCluster.size();
           int b = sampleCluster.size() - a;
           // Note: all hypermutated samples are resequenced, and in our clusters.
           int c = hyperMutatedSamples.size() - a;
           int d = samples.length - sampleCluster.size() - c;
           double pvalue = fisher.getTwoTailedP(a,
                                                b,
                                                c,
                                                d);
           System.out.printf("Contingence table: %d, %d, %d, %d%n", a, b, c, d);
           System.out.println("pvalue from Fisher: " + pvalue);
           pvalue = MathUtilities.calculateHypergeometricPValue(samples.length, 
                                                                sampleCluster.size(), 
                                                                hyperMutatedSamples.size(), 
                                                                hyperInCluster.size());
           if (pvalue > 0.50) // Try to see another way
               pvalue = MathUtilities.calculateHypergeometricPValue(samples.length, 
                                                                    sampleCluster.size(), 
                                                                    samples.length - hyperMutatedSamples.size(), 
                                                                    sampleCluster.size() - hyperInCluster.size());
           System.out.println("pvalue from hyper: " + pvalue);
       }       
       // Check sample survival rate
       Map<String, Integer> sampleToSurvival = loadSampleToSurvivalRate();
       // Get the average survival rate
       int total = 0;
       for (String sample : sampleToSurvival.keySet())
           total += sampleToSurvival.get(sample);
       double avgRate = (double) total / sampleToSurvival.size();
       System.out.println("\nAverage survival rate: " + avgRate);
       index = 1;
       // Check survival rate for each cluster
       DescriptiveStatistics stat = new DescriptiveStatistics();
       for (List<String> sampleCluster : sampleClusters) {
           System.out.printf("%nCluster %d: %d%n",
                             index ++,
                             sampleCluster.size());
           stat.clear();
           // Want to do a t-test
           double[] values = new double[sampleCluster.size()];
           for (int i = 0; i < sampleCluster.size(); i++) {
               values[i] = sampleToSurvival.get(sampleCluster.get(i));
               stat.addValue(values[i]);
           }
           System.out.println("Average in the cluster: " + stat.getMean());
           double ttest = TestUtils.tTest(avgRate, values);
           System.out.println("pvalue: " + ttest);
       }
//       // Want to generate a text file for loading into R
//       System.out.println("Sample\tVital\tLength\thCluster\tRecurrence");
//       index = 1;
//       Map<String, Integer> sampleToVital = loadSampleToVitalStatus();
//       for (List<String> sampleCluster : sampleClusters) {
//           for (String sample : sampleCluster) {
//               int vital = sampleToVital.get(sample);
//               int length = sampleToSurvival.get(sample);
//               String recurrence = sampleToRec.get(sample);
//               System.out.println(sample + "\t" + 
//                                  vital + "\t" +
//                                  length + "\t" + 
//                                  index + "\t" + 
//                                  recurrence);
//           }
//           index ++;
//       }
    }
    
    @Test
    public void checkSamplesInNetworkClusters() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
//        // The following is used to pick up CNV genes for resequenced samples
//        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToCNVGenes(samples);
//        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
//        Map<String, Integer> sampleToSurvival = loadSampleToTimeSpan();
        // Check for science altered genes
        //Map<String, Set<String>> sampleToAlteredGenes = loadAllScienceSampleToMutatedGenes();
//        Map<String, Set<String>> sampleToMutatedGenes = loadAllScienceSampleToMutatedGenes();
//        Map<String, Set<String>> sampleToCNVGenes = loadAllScienceSampleToCNVGenes();
//        Map<String, Set<String>> sampleToAlteredGenes = new HashMap<String, Set<String>>();
//        // Use the following statements to get samples having both CNV and mutation data
//        // since CNV data is available to the 22 discovery screen samples only. The following
//        // statements are used for filtering mutated genes only
//        for (String sample : sampleToCNVGenes.keySet()) {
//            Set<String> genes = sampleToMutatedGenes.get(sample);
//            if (genes == null)
//                continue;
//            sampleToAlteredGenes.put(sample, genes);
//        }
        // Altered genes for 22 discovery genes only.
//        Map<String, Set<String>> sampleToAlteredGenes = getScienceGBMSampleToAlteredGenes();

        System.out.println("Total samples: " + sampleToAlteredGenes.size());
        //String clusterFileName = TCGA_GBM_DIR + "ClustersFIsInMutatedGenesFromTwoGBMs071509_010.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInFIsInHomoCNVMutatedGenes070709.txt";
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesSpectral032910.txt";
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes_2Samples_093009.txt";
//        String clusterFileName = TCGA_GBM_DIR + "ClustersInScienceGBMAlteredGenes092909.txt";
     //   String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenesInPathwayFIs080509_115.txt";
//        List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(clusterFileName);
        
        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToAlteredGenes);
        NetworkModularityCalculator calculator = new NetworkModularityCalculator();
        List<NetworkModule> modules = calculator.calculateClusterModularityWithFDR(allGenes,
                                                                                   true,
                                                                                   0);
        List<Set<String>> clusters = calculator.convertModuleObjectToSet(modules);
        System.out.println("Module\tSize");
        for (int i = 0; i < clusters.size(); i++)
            System.out.println(i + "\t" + clusters.get(i).size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        clusterAnalyzer.annotateNetworkClusters(clusters, 10);
        
        System.out.println("Total clusters: " + clusters.size());
        // Check how many samples having 0 and 1 clusters
        int counter = 0;
        System.out.println("Sample\tClusters");
        int cluster0 = 0;
        int cluster1 = 0;
        List<String> posSamples = new ArrayList<String>();
        List<String> negSamples = new ArrayList<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> genes = sampleToAlteredGenes.get(sample);
            // Want to get the list of clusters
            List<Integer> clusterIds = new ArrayList<Integer>();
            List<String> sharedGenes = new ArrayList<String>();
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> set = clusters.get(i);
                for (String gene : genes) {
                    if (set.contains(gene)) {
                        clusterIds.add(i);
                        break;
                    }
                }
                Set<String> shared = new HashSet<String>(set);
                shared.retainAll(genes);
                sharedGenes.add(StringUtils.join(":", new ArrayList<String>(shared)));
            }
            //System.out.println(sample + "\t" + clusterIds + "\t" + sharedGenes);
            System.out.println(sample + "\t" + clusterIds);
            if (isSampleInClusters(clusterIds)) {
                counter ++;
                posSamples.add(sample);
            }
            else
                negSamples.add(sample);
            if (clusterIds.contains(0))
                cluster0 ++;
            if (clusterIds.contains(2))
                cluster1 ++;
        }
        double percent = (double) counter / sampleToAlteredGenes.size();
        System.out.println("Total samples checked: " + sampleToAlteredGenes.size());
        System.out.println("Total samples having 0 and 2: " + counter + "(" + percent + ")");
        System.out.println("Total samples in cluster 0: " + cluster0);
        System.out.println("Total samples in cluster 2: " + cluster1);
//        // Print out samples
//        Collections.sort(posSamples);
//        Collections.sort(negSamples);
//        System.out.println("\nSamples have altered genes in both clusters 0 and 1: " + posSamples.size());
//        for (String sample : posSamples)
//            System.out.println(sample);
//        System.out.println("\nSamples have altered genes not in both clusters 0 and 1: " + negSamples.size());
//        for (String sample : negSamples)
//            System.out.println(sample);
//        // Do a recurrence enrichment analysis
//        Map<String, String> sampleToRec = loadSampleToRecurrent();
//        Set<String> noSamplesInPos = new HashSet<String>();
//        for (String sample : posSamples) {
//            String rec = sampleToRec.get(sample);
//            if (rec.equals("No"))
//                noSamplesInPos.add(sample);
//        }
//        int a = noSamplesInPos.size();
//        int b = posSamples.size() - a;
//        // Total rec or sec: 19
//        int c = sampleToAlteredGenes.size() - 19 - a;
//        int d = negSamples.size() - c;
//        System.out.printf("%nContigence table: %d, %d, %d, %d%n",
//                          a, b, c, d);
//        double pvalue = new FisherExact(100).getTwoTailedP(a, b, c, d);
//        System.out.println("pvalue from Fisher: " + pvalue);
//        pvalue = MathUtilities.calculateHypergeometricPValue(sampleToAlteredGenes.size(), 
//                                                             posSamples.size(),
//                                                             sampleToAlteredGenes.size() - 19,
//                                                             noSamplesInPos.size());
//        System.out.println("pvalue from hyper: " + pvalue);
    }

    private Map<String, Set<String>> loadSampleToCNVGenes(List<String> samples)
            throws Exception {
        Map<String, Set<String>> sampleToAlteredGenes = loadSampleToCNVGenes();
        for (Iterator<String> it = sampleToAlteredGenes.keySet().iterator(); it.hasNext();) {
            String sample = it.next();
            if (samples.contains(sample))
                continue;
            it.remove();
        }
        return sampleToAlteredGenes;
    }
    
    @Test
    public void generateFIsAmongMutatedGenesAndCNVs() throws Exception {
        List<String> samples = loadResequencedSamples();
        Set<String> genes = getAlteredGenesInSamples(samples);
        System.out.println("Total altered genes: " + genes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        genes.retainAll(totalGenes);
        System.out.println("Genes in FI network: " + genes.size());
        Set<String> neededFIs = InteractionUtilities.getFIs(genes, fis);
        Set<String> tmp = InteractionUtilities.grepIDsFromInteractions(neededFIs);
        System.out.println("Genes after FI filtering: " + tmp.size());
        //String outFileName = TCGA_GBM_DIR + "FIsInHomoCNVMutatedGenes.txt";
        //fu.saveInteractions(neededFIs, outFileName);
    }

    Set<String> getAlteredGenesInSamples(List<String> samples) throws Exception {
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        Set<String> genes = new HashSet<String>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> set = sampleToAlteredGenes.get(sample);
            genes.addAll(set);
        }
        return genes;
    }
    
    @Test
    public void checkCNVs() throws IOException {
        String name = "CNV_Amplication_RAE.txt";
        Set<String> raeAmpGenes = loadCNVGenes(TCGA_GBM_DIR + name);
        name = "CNV_Deletion_RAE.txt";
        Set<String> raeDeleteGenes = loadCNVGenes(TCGA_GBM_DIR + name);
        String[] ampFileNames = new String[] {
                "CNV_Amplication_GISTC.txt",
                "CNV_Amplication_RAE.txt",
                "CNV_Amplication_GTS.txt"
        };
        String[] deleteFileNames = new String[] {
                "CNV_Deletion_GISTC.txt",
                "CNV_Deletion_RAE.txt",
                "CNV_Deletion_GTS.txt"
        };
        for (String fileName : ampFileNames) {
            Set<String> genes = loadCNVGenes(TCGA_GBM_DIR + fileName);
            System.out.println(fileName + ": " + genes.size());
            genes.removeAll(raeAmpGenes);
            System.out.println("Not shared with RAE: " + genes.size());
        }
        for (String fileName : deleteFileNames) {
            Set<String> genes = loadCNVGenes(TCGA_GBM_DIR + fileName);
            System.out.println(fileName + ": " + genes.size());
            genes.removeAll(raeDeleteGenes);
            System.out.println("Not shared with RAE: " + genes.size());
        }
    }
    
    private Set<String> loadCNVGenes(String fileName) throws IOException {
        fu.setInput(fileName);
        Set<String> geneSet = new HashSet<String>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            String[] tokens = line.split(",");
            for (String token : tokens)
                geneSet.add(token.trim());
        }
        fu.close();
        return geneSet;
    }
    
    public Map<String, Set<String>> loadSampleToMutations() throws Exception {
        String fileName = TCGA_GBM_DIR + "GlioblastomaMutationTable.txt";
        return loadSampleToMutatedGenes(fileName);
    }

    protected Map<String, Set<String>> getSampleToAlteredGenes(List<String> samples) throws Exception {
        Map<String, Set<String>> sampleToCNVGenes = loadSampleToCNVGenes();
        Map<String, Set<String>> sampleToMutations = loadSampleToMutations();
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        for (String sample : samples) {
            Set<String> cnvGenes = sampleToCNVGenes.get(sample);
            if (cnvGenes == null)
                cnvGenes = new HashSet<String>();
            Set<String> mutationGenes = sampleToMutations.get(sample);
            if (mutationGenes == null)
                mutationGenes = new HashSet<String>();
            Set<String> all = new HashSet<String>();
            all.addAll(cnvGenes);
            all.addAll(mutationGenes);
            sampleToGenes.put(sample, all);
        }
        return sampleToGenes;
    }
    
    private Map<String, Double> calculateSamplePairToScoreFromShare(List<String> samples,
                                                                    Map<String, Set<String>> sampleToAlteredGenes) throws Exception {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Double> pairToScore = new HashMap<String, Double>();
        for (int i = 0; i < samples.size() - 1; i++) {
            String sample1 = samples.get(i);
            Set<String> genes1 = sampleToAlteredGenes.get(sample1);
            genes1.retainAll(allGenes);
            Set<String> fis1 = getFIs(genes1, fis);
            for (int j = i + 1; j < samples.size(); j++) {
                String sample2 = samples.get(j);
                Set<String> genes2 = sampleToAlteredGenes.get(sample2);
                genes2.retainAll(allGenes);
                Set<String> fis2 = getFIs(genes2, fis);
                // Calculate score
                // Score for shared genes
//                int geneShared = sharedNumber(genes1, genes2);
//                double genePvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(), 
//                                                                            genes1.size(),
//                                                                            genes2.size(),
//                                                                            geneShared);
                int fiShared = sharedNumber(fis1, fis2);
//                double fiPValue = MathUtilities.calculateHypergeometricPValue(fis.size(), 
//                                                                              fis1.size(),
//                                                                              fis2.size(), 
//                                                                              fiShared);
                // Use binomail for fast approximate calculation
                double ratio = (double) fis1.size() / fis.size();
                double fiPValue = MathUtilities.calculateBinomialPValue(ratio, fis2.size(), fiShared);
                pairToScore.put(sample1 + "\t" + sample2, fiPValue);
                System.out.println(sample1 + "\t" + sample2 + ": " + fiPValue);
            }
        }
        return pairToScore;
    }
    
    private Set<String> getFIs(Set<String> genes,
                               Set<String> fis) {
        Set<String> rtn = new HashSet<String>();
        for (String gene : genes) {
            Set<String> tmp = InteractionUtilities.grepFIsContains(gene, fis);
            rtn.addAll(tmp);
        }
        return rtn;
    }
                
    private int sharedNumber(Set<String> set1, Set<String> set2) {
        int c = 0;
        for (String gene1 : set1) {
            if (set2.contains(gene1))
                c ++;
        }
        return c;
    }
    
    
     
    /**
     * Cluster GBM samples based on pair-wise shortest distances among altered genes in samples.
     * @throws Exception
     */
    @Test
    public void hierachicalClusteringSamples() throws Exception {
        List<String> resequencedSamples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(resequencedSamples);
        //Map<String, Set<String>> sampleToAlteredGenes = loadSampleToMutations();
        resequencedSamples.retainAll(sampleToAlteredGenes.keySet());
        hierchicallyClusterSamples(resequencedSamples, sampleToAlteredGenes);
    }
    
    /**
     * This method is used to generate a matrix from samples to pathways based on connection
     * degrees in pathways using FIs.
     * @throws Exception
     */
    @Test
    public void generateSampleToNetworkModulesInDegrees() throws Exception {
        //Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        // For TCGA samples
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String outFileName = TCGA_GBM_DIR + "TCGASampleToMCLNetworkModulesConnDegreeNorm120709.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        String clusterFileName = TCGA_GBM_DIR + "MCLForFIsInTCGAAlteredGenes071609_120709.txt";
        clusterAnalyzer.generateSampleToModulesMatrixInDegrees(sampleToAlteredGenes, 
                                                               clusterFileName, 
                                                               outFileName,
                                                               2); // Pick modules at least 2 or more    
    }
    
    @Test
    public void generateSampleToNetworkModuleMatrix() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String outFileName = TCGA_GBM_DIR + "TCGASampleToMCLNetworkModulesBinary120709.txt";
        //String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        String clusterFileName = TCGA_GBM_DIR + "MCLForFIsInTCGAAlteredGenes071609_120709.txt";
        clusterAnalyzer.generateSampleToModulesMatrix(sampleToAlteredGenes, 
                                                      clusterFileName, 
                                                      outFileName,
                                                      2);
    }

    /**
     * Cluster GBM samples based on alteration occurences in the network clusters.
     * @throws Exception
     */
    @Test
    public void hierarchicalClusterSamplesOnNetworkClusters() throws Exception {
        String clusterFileName = TCGA_GBM_DIR + "ClustersInTCGAAlteredGenes071609.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        // For TCGA samples
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        // Want to create a sample to a boolean vector map
        Map<String, List<Boolean>> sampleToVector = new HashMap<String, List<Boolean>>();
        for (String sample : sampleToAlteredGenes.keySet()) {
            Set<String> alteredGenes = sampleToAlteredGenes.get(sample);
            List<Boolean> vector = new ArrayList<Boolean>();
            for (Set<String> cluster : clusters) {
                //if (cluster.size() <= 5)
                //    continue;
                Set<String> shared = InteractionUtilities.getShared(cluster, alteredGenes);
                if (shared.size() > 0)
                    vector.add(Boolean.TRUE);
                else
                    vector.add(Boolean.FALSE);
                //vector.add(shared.size());
            }
            System.out.println(sample + ": " + StringUtils.join(", ", vector));
            sampleToVector.put(sample, vector);
        }
        // Create an output
        //String outFileName = TCGA_GBM_DIR + "TCGASampleToNetworkClusters.txt";
        // Have to place the following statement before the output statements,
        // otherwise, it will not work.
        Map<String, String> sampleToRecurrence = loadSampleToRecurrent();
        String outFileName = TCGA_GBM_DIR + "TCGASampleToNetworkClustersNoSamples.txt";
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Total cluster
        int totalCluster = clusters.size();
        builder.append("Sample");
        for (int i = 0; i < totalCluster; i++)
            builder.append("\tModule").append(i);
        fu.printLine(builder.toString());
        for (String sample : sampleToVector.keySet()) {
            String recurrence = sampleToRecurrence.get(sample);
            if (!recurrence.equals("No"))
                continue;
            builder.setLength(0);
            builder.append(sample);
            List<Boolean> vector = sampleToVector.get(sample);
            for (Boolean b : vector) {
                builder.append("\t");
                if (b)
                    builder.append(1);
                else
                    builder.append(0);
            }
            //System.out.println(builder.toString());
            fu.printLine(builder.toString());
        }
        fu.close();
        if (true)
            return;
        // Create pair-wise distances among samples
        Map<String, Double> samplePairToDist = new HashMap<String, Double>();
        List<String> sampleList = new ArrayList<String>(sampleToAlteredGenes.keySet());
        for (int i = 0; i < sampleList.size() - 1; i++) {
            String sample1 = sampleList.get(i);
            List<Boolean> vector1 = sampleToVector.get(sample1);
            for (int j = i + 1; j < sampleList.size(); j++) {
                String sample2 = sampleList.get(j);
                List<Boolean> vector2 = sampleToVector.get(sample2);
                //double dist = MathUtilities.calculateHammingDistance(vector1, vector2);
                double dist = MathUtilities.calculateNetworkDistance(vector1,
                                                                     vector2, 
                                                                     clusters);
                samplePairToDist.put(sample1 + "\t" + sample2,
                                     dist);
            }
        }
        // Do a hierarchical clustering
        List<HierarchicalClusterNode> sampleClusters = hierarchicalCluster(sampleList, 
                                                               samplePairToDist,
                                                               false);
        HierarchicalClusterNode root = sampleClusters.get(0);
        outputCluster(root);
    }
    
    protected double calculateAverageDistance(Set<String> set1,
                                              Set<String> set2,
                                              BreadthFirstSearch bfs,
                                              Map<String, Set<String>> idToPartners) {
        List<String> list1 = new ArrayList<String>(set1);
        List<String> list2 = new ArrayList<String>(set2);
        int total = 0;
        int count = 0;
        for (String gene1 : list1) {
            Map<String, Integer> geneToDistances = bfs.getDistances(gene1, 
                                                                    list2,
                                                                    idToPartners);
            count += geneToDistances.size();
            for (Integer p : geneToDistances.values())
                total += p;
        }
        return (double) total / count;
    }
    
    private double calculateMinShortestPathOneWay(Set<String> source,
                                                  Set<String> target,
                                                  BreadthFirstSearch bfs,
                                                  Map<String, Set<String>> idToPartners) {
        List<String> list1 = new ArrayList<String>(source);
        List<String> list2 = new ArrayList<String>(target);
        int total = 0;
        int count = 0;
        // From list1 to list2
        for (String gene1 : list1) {
            Map<String, Integer> geneToDistances = bfs.getDistances(gene1, 
                                                                    list2,
                                                                    idToPartners);
            // Find the shortest path
            int shortest = Integer.MAX_VALUE;
            for (Integer tmp : geneToDistances.values()) {
                if (tmp < shortest)
                    shortest = tmp;
            }
            total += shortest;
            count ++;
        }
        return (double) total / count;
    }
    
    @Test
    public void checkShortestPathInSamples() throws Exception {
        Map<String, Set<String>> sampleToCNVGenes = loadSampleToCNVGenes();
        System.out.println("Total samples: " + sampleToCNVGenes.size());
        List<String> resequencedSamples = loadResequencedSamples();
        System.out.println("total resequenced: " + resequencedSamples.size());
        Map<String, Set<String>> sampleToMutations = loadSampleToMutations();
        System.out.println("Total mutation samples: " + sampleToMutations.size());
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        System.out.println();
        System.out.println("Index\tSample\tCNV_Genes\tSize\tShortestPath\t" +
        		           "Mutation_Genes\tSize\tShortestPath\t" +
        		           "CNV_Mutation_Genes\tSize\tShortestPath");
        int c = 1;
        for (String sample : resequencedSamples) {
            Set<String> cnvGenes = sampleToCNVGenes.get(sample);
            if (cnvGenes == null)
                cnvGenes = new HashSet<String>();
            double cnvPath = calculateAverageShortestPath(bfs, allGenes, nodeToEdges, cnvGenes);
            Set<String> mutations = sampleToMutations.get(sample);
            if (mutations == null)
                mutations = new HashSet<String>();
            double mutationPath = calculateAverageShortestPath(bfs, allGenes, nodeToEdges, mutations);
            Set<String> all = new HashSet<String>(cnvGenes);
            all.addAll(mutations);
            double allPath = calculateAverageShortestPath(bfs, allGenes, nodeToEdges, all);
            System.out.println(c + "\t" +
                               sample + "\t" + 
                               cnvGenes + "\t" +
                               cnvGenes.size() + "\t" +
                               cnvPath + "\t" + 
                               mutations + "\t" +
                               mutations.size() + "\t" + 
                               mutationPath + "\t" +
                               all + "\t" +
                               all.size() + "\t" + 
                               allPath);
            c ++;
        }
    }

    private double calculateAverageShortestPath(BreadthFirstSearch bfs,
                                                Set<String> allGenes,
                                                Map<TreeNode, List<Edge>> nodeToEdges,
                                                Set<String> cnvGenes) {
        double cnvPath = -1.0;
        List<String> list = new ArrayList<String>(cnvGenes);
        list.retainAll(allGenes);
        if (list.size() > 2) {
            cnvPath = calculateShortestPath(list, bfs, nodeToEdges);
        }
        return cnvPath;
    }
    
    @Test
    public void testLoadSampleToMutations() throws Exception {
        Map<String, Set<String>> sampleToCNVGenes = loadSampleToCNVGenes();
        System.out.println("Total samples: " + sampleToCNVGenes.size());
        List<String> resequencedSamples = loadResequencedSamples();
        System.out.println("total resequenced: " + resequencedSamples.size());
        Map<String, Set<String>> sampleToMutations = loadSampleToMutations();
        System.out.println("Total mutation samples: " + sampleToMutations.size());
        System.out.println();
        System.out.println("Index\tSample\tCNV_Genes\tMutation_Genes\tCNV_Mutation_Genes");
        int c = 1;
        Set<String> grandAll = new HashSet<String>();
        for (String sample : resequencedSamples) {
            Set<String> cnvGenes = sampleToCNVGenes.get(sample);
            if (cnvGenes == null)
                cnvGenes = new HashSet<String>();
            Set<String> mutations = sampleToMutations.get(sample);
            if (mutations == null)
                mutations = new HashSet<String>();
            Set<String> all = new HashSet<String>(cnvGenes);
            all.addAll(mutations);
            System.out.println(c + "\t" +
                               sample + "\t" + 
                               cnvGenes.size() + "\t" +
                               mutations.size() + "\t" + 
                               all.size());
            c ++;
            grandAll.addAll(all);
        }
        System.out.println("\nAll CNV and mutated genes: " + grandAll.size());
    }
    
    public Map<String, Set<String>> loadSampleToCNVGenes() throws Exception {
        Map<String, Set<CopyNumberVariation>> sampleToCNVs = loadSampleToCNVs();
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        for (String sample : sampleToCNVs.keySet()) {
            Set<CopyNumberVariation> cnvs = sampleToCNVs.get(sample);
            Set<String> genes = new HashSet<String>();
            for (CopyNumberVariation cnv : cnvs) {
                genes.add(cnv.getGene());
            }
            // Use sample ID only
            sampleToGenes.put(sample.substring(0, 12), genes);
        }
        return sampleToGenes;
    }
    
    public Map<String, Set<String>> loadSampleToCNVAmplifiedGenes() throws Exception {
        return loadSampleToCNVGenes(true);
    }
    
    public Map<String, Set<String>> loadSampleToCNVDeletedGenes() throws Exception {
        return loadSampleToCNVGenes(false);
    }

    private Map<String, Set<String>> loadSampleToCNVGenes(boolean isForAmplied) throws Exception {
        Map<String, Set<CopyNumberVariation>> sampleToCNVs = loadSampleToCNVs();
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        for (String sample : sampleToCNVs.keySet()) {
            Set<CopyNumberVariation> cnvs = sampleToCNVs.get(sample);
            Set<String> genes = new HashSet<String>();
            for (CopyNumberVariation cnv : cnvs) {
                // For amplified genes
                if (isForAmplied && cnv.getValue() == 2)
                    genes.add(cnv.getGene());
                // For deleted genes
                else if (!isForAmplied && cnv.getValue() == -2)
                    genes.add(cnv.getGene());
            }
            // Use sample ID only
            sampleToGenes.put(sample.substring(0, 12), genes);
        }
        return sampleToGenes;
    }
    
    
    
    public Map<String, Set<CopyNumberVariation>> loadSampleToCNVs() throws Exception {
        // There is an error in the original file: TCGA-GBM-RAE-genemap-n216-20080510-dscrt.txt.
        // Change the CB column sample barcode from TCGA-06-0178-10B-01D to TCGA-06-0178-01A-01D.
        // 10B should be code for normal sample, which is not correct in this file.
        String fileName = TCGA_GBM_DIR + "TCGA-GBM-RAE-genemap-n216-20080510-dscrt_1.txt";
        fu.setInput(fileName);
        Map<String, Set<CopyNumberVariation>> sampleToCNVs = new HashMap<String, Set<CopyNumberVariation>>();
        String line = fu.readLine();
        // Get the samples
        List<String> samples = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++) {
            if (tokens[i].equals("Chr"))
                break;
            samples.add(tokens[i]);
        }
        int index = 0;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            index = tokens[0].indexOf(":");
            String gene = tokens[0].substring(index + 1);
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("0") ||
                    tokens[i].equals("NA"))
                    continue; // Nothing useful here
                if (i > samples.size())
                    break;
                // Remember to deduct 1 in the index
                String sample = samples.get(i - 1);
                int value = Integer.parseInt(tokens[i]);
                CopyNumberVariation cnv = new CopyNumberVariation();
                cnv.setGene(gene);
                cnv.setValue(value);
                Set<CopyNumberVariation> cnvs = sampleToCNVs.get(sample);
                if (cnvs == null) {
                    cnvs = new HashSet<CopyNumberVariation>();
                    sampleToCNVs.put(sample, cnvs);
                }
                cnvs.add(cnv);
            }
        }
        fu.close();
        // Do not use the following expression validation since the another two filters,
        // CNV_Amplifcation_RAE and CNV_Deletion_RAE have considered gene expression already.
        //Map<String, List<GeneExpressionValue>> sampleToGeneExpValues = loadSampleToGeneExpValues();
        // Want to clean up CNVs a little: if a gene is deleted and amplified in the same sample
        // don't use it.
        Set<String> cnvAmplificatedGenes = loadCNVGenes(TCGA_GBM_DIR + "CNV_Amplication_RAE.txt");
        //System.out.println("Total amplified genes: " + cnvAmplificatedGenes.size());
        Set<String> cnvDeletedGenes = loadCNVGenes(TCGA_GBM_DIR + "CNV_Deletion_RAE.txt");
        //System.out.println("Total deleted genes: " + cnvDeletedGenes.size());
        for (Iterator<String> it = sampleToCNVs.keySet().iterator(); it.hasNext();) {
            String sample = it.next();
            Set<CopyNumberVariation> cnvs = sampleToCNVs.get(sample);
//            List<GeneExpressionValue> expValues = getMatchedGeneExpressionValues(sample, 
//                                                                                 sampleToGeneExpValues);
//            if (expValues == null)
//                System.out.println(sample + " has no expression values!");
            cleanUpCNVs(cnvs,
                        null,
                        cnvAmplificatedGenes, 
                        cnvDeletedGenes);
            //if (cnvs.size() == 0)
            //    it.remove();
        }
        return sampleToCNVs;
    }
    
    /**
     * Sample from CNV may have different barcode from sample from gene expression. So a mapping should
     * be done here to get gene expression for the same cancer sample.
     * @param sample
     * @param sampleToGeneExpValues
     * @return
     */
    private List<GeneExpressionValue> getMatchedGeneExpressionValues(String sample,
                                                                     Map<String, List<GeneExpressionValue>> sampleToGeneExpValues) {
        // A barcode is something like this: TCGA-02-0001-01C-01D. The most important one is: TCGA-02-0001-01
        String sampleId = sample.substring(0, sample.length() - 5);
        for (String tmp : sampleToGeneExpValues.keySet()) {
            if (tmp.startsWith(sampleId))
                return sampleToGeneExpValues.get(tmp);
        }
        return null;
    }
    
    private Map<String, List<GeneExpressionValue>> loadSampleToGeneExpValues() throws Exception {
        long time1 = System.currentTimeMillis();
        String fileName = TCGA_GBM_DIR + "GeneExp_z_scores_051409.txt";
        fu.setInput(fileName);
        List<String> controls = getControlSamples();
        List<String> cancerSamples = new ArrayList<String>();
        List<String> controlSamples = new ArrayList<String>();
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        //System.out.println("Total samples: " + (tokens.length - 1));
        List<String> all = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++) {
            String sample = tokens[i];
            sample = sample.replaceAll("\\.", "-");
            all.add(sample);
            if (isControlSample(sample, controls)) 
                controlSamples.add(sample);
            else 
                cancerSamples.add(sample);
        }
        //System.out.println("Total cancer samples: " + samples.size());
        Map<String, List<GeneExpressionValue>> sampleToGeneExpValues = new HashMap<String, List<GeneExpressionValue>>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Want to get control values to calculate t-stat
            List<Double> controlValues = new ArrayList<Double>();
            // The following loop is used to get all control values
            for (int i = 1; i < tokens.length; i++) {
                String sample = all.get(i - 1);
                if (controlSamples.contains(sample)) {
                    controlValues.add(new Double(tokens[i]));
                }
            }
            // Required by the math lib
            double[] controlValues1 = new double[controlValues.size()];
            for (int i = 0; i < controlValues.size(); i++)
                controlValues1[i] = controlValues.get(i);
            // The following loop is used to calculate expression values
            for (int i = 1; i < tokens.length; i++) {
                String sample = all.get(i - 1);
                if (cancerSamples.contains(sample)) {
                    List<GeneExpressionValue> values = sampleToGeneExpValues.get(sample);
                    if (values == null) {
                        values = new ArrayList<GeneExpressionValue>();
                        sampleToGeneExpValues.put(sample, values);
                    }
                    GeneExpressionValue value = new GeneExpressionValue();
                    value.setGene(tokens[0]);
                    // Use in a reverse way
                    double zscore = Double.parseDouble(tokens[i]);
                    value.setTvalue(-TestUtils.t(zscore, controlValues1));
                    value.setPvalue(TestUtils.tTest(zscore, controlValues1));
                    value.setZvalue(zscore);
                    values.add(value);
                }
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for loading gene expression values: " + (time2 - time1));
        return sampleToGeneExpValues;
    }
    
    private void cleanUpCNVs(Set<CopyNumberVariation> cnvs,
                             List<GeneExpressionValue> expValues,
                             Set<String> cnvAmplifiedGenes,
                             Set<String> cnvDeleteGenes) {
        Map<String, Set<String>> geneToTypes = new HashMap<String, Set<String>>();
        for (CopyNumberVariation cnv : cnvs) {
            String gene = cnv.getGene();
            Set<String> types = geneToTypes.get(gene);
            if (types == null) {
                types = new HashSet<String>();
                geneToTypes.put(gene, types);
            }
            int value = cnv.getValue();
            if (value > 0)
                types.add("+");
            else if (value < 0)
                types.add("-");
        }
        // Clean up: keep only amplified or deleted, but not both.
        for (Iterator<CopyNumberVariation> it = cnvs.iterator(); it.hasNext();) {
            CopyNumberVariation cnv = it.next();
            Set<String> types = geneToTypes.get(cnv.getGene());
            if (types.size() > 1)
                it.remove();
        }
        // Clean up based on gene exp value
        // Use gene expression cannot help a lot
        if (expValues != null) {
            double pvalue_cutoff = 0.0001;
            for (Iterator<CopyNumberVariation> it = cnvs.iterator(); it.hasNext();) {
                CopyNumberVariation cnv = it.next();
                // Check its gene expression
                GeneExpressionValue geneExp = null;
                for (GeneExpressionValue value : expValues) {
                    if (value.getGene().equals(cnv.getGene())) {
                        geneExp = value;
                        break;
                    }
                }
                if (geneExp == null)
                    it.remove();
                else {
                    // One up one down
                    if (cnv.getValue() * geneExp.getTvalue() < 0) 
                        it.remove();
                    else if (geneExp.getPvalue() > pvalue_cutoff)
                        it.remove();
                }
            }
        }
        // Do another filter based on a global analysis
        for (Iterator<CopyNumberVariation> it = cnvs.iterator(); it.hasNext();) {
            CopyNumberVariation cnv = it.next();
            if (cnv.getValue() > 0 && cnvAmplifiedGenes.contains(cnv.getGene()))
                continue;
            if (cnv.getValue() < 0 && cnvDeleteGenes.contains(cnv.getGene()))
                continue;
            it.remove();
        }
        // Use homo-deletions and multiple-amplication CNVs only
        for (Iterator<CopyNumberVariation> it = cnvs.iterator(); it.hasNext();) {
            CopyNumberVariation cnv = it.next();
            if (cnv.getValue() == 2 || cnv.getValue() == -2)
                continue;
            it.remove();
        }
    }
    
    List<String> loadResequencedSamples() throws IOException {
        String fileName = TCGA_GBM_DIR + "IndividualSamples.txt";
        fu.setInput(fileName);
        List<String> samples = new ArrayList<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            //System.out.println(line);
            String[] tokens = line.split("\t");
            if (tokens.length < 9)
                continue;
            if (tokens[8].equals("Yes")) {
                //if (tokens.length > 9 && tokens[9].equals("Yes"))
                //    continue;
                samples.add(tokens[0]);
            }
        }
        fu.close();
        return samples;
    }
    
    /**
     * Used to check linker genes.
     */
    @Test
    public void checkGBMLinkerGenes() throws Exception {
        Set<String> fis = fu.loadInteractions(DIR_NAME + "FIsInNautreGBM050409.txt");
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        fiGenes.removeAll(mutatedGenes);
        List<String> linkerGenes = new ArrayList<String>(fiGenes);
        // Sorting based on gene rankers from MSKCC
        final Map<String, Double> geneToScore = loadGeneRankers();
        Collections.sort(linkerGenes, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double score1 = geneToScore.get(gene1);
                if (score1 == null)
                    score1 = 0.0d;
                Double score2 = geneToScore.get(gene2);
                if (score2 == null)
                    score2 = 0.0d;
                return score2.compareTo(score1);
            }
        });
        System.out.println("Total linker genes: " + linkerGenes.size());
        // Want to print out protein names
        Map<String, String> geneNameToProteinName = new HibernateFIReader().loadShortNameToNameMap(linkerGenes);
        Map<String, Double> geneToTValue = new HashMap<String, Double>();
        Map<String, Double> geneToPValue = new HashMap<String, Double>();
        String fileName = TCGA_GBM_DIR + "GeneExp_t_value.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToTValue.put(tokens[0],
                             new Double(tokens[1])); 
            geneToPValue.put(tokens[0],
                             new Double(tokens[2]));
        }
        System.out.println("Gene_Name\tRanker_Score\tt-value\tp-value\tCNV_Amplificated\tCNV_Deleted\tProtein_Name");
        Set<String> cnvAmplificatedGenes = loadCNVGenes(TCGA_GBM_DIR + "CNV_Amplication.txt");
        System.out.println("Total amplicated genes: " + cnvAmplificatedGenes.size());
//        for (String gene : cnvAmplificatedGenes)
//            System.out.println(gene);
        Set<String> cnvDeletedGenes = loadCNVGenes(TCGA_GBM_DIR + "CNV_Deletion.txt");
        System.out.println("Total deleted genes: " + cnvDeletedGenes.size());
        for (String linker : linkerGenes) {
            Double score = geneToScore.get(linker);
            String proteinName = geneNameToProteinName.get(linker);
            if (proteinName == null)
                proteinName = "";
            String tvalue = null;
            if (geneToTValue.get(linker) == null)
                tvalue = "";
            else
                tvalue = String.format("%.2e", geneToTValue.get(linker));
            String pvalue = null;
            if (geneToPValue.get(linker) == null)
                pvalue = "";
            else
                pvalue = String.format("%.2e", geneToPValue.get(linker));
            boolean isCNVAmplificated = cnvAmplificatedGenes.contains(linker);
            boolean isCNVDeleted = cnvDeletedGenes.contains(linker);
            System.out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                              linker,
                              (score == null ? "" : score.toString()),
                              tvalue,
                              pvalue, 
                              isCNVAmplificated + "",
                              isCNVDeleted + "",
                              proteinName);
        }
    }
    
    @Test
    public void fixGeneRankerFile() throws IOException {
        String inFileName = TCGA_GBM_DIR + "Gene_Ranker_All.txt";
        String outFileName = TCGA_GBM_DIR + "Gene_Ranker_All_Fixed.txt";
        fu.setInput(inFileName);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        String line = fu.readLine();
        outFu.printLine(line);
        String preLine = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Gene")) {
                outFu.printLine(preLine + "\t" + line);
            }
            else if (line.startsWith("Num"))
                continue;
            else
                preLine = line;
        }
        fu.close();
        outFu.close();
    }
    
    /**
     * This method is used to check the degree of genes connecting to the mutated genes.
     * @throws IOException
     */
    @Test
    public void generateConnectionToMutatedGeneDegreeFile() throws IOException {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> fiToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        Set<String> allGenes = new HashSet<String>(fiToPartners.keySet());
        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        String outFileName = TCGA_GBM_DIR + "GeneToFIMutatedDegree.txt";
        fu.setOutput(outFileName);
        fu.printLine("Gene\tMutatedConnectionDegree\ttvalue");
        String geneExpFileName = "datasets/TCGA/GBM/GeneExp_t_value.txt";
        Map<String, Double> geneToExpTValue = loadGeneExpTValue(geneExpFileName);
        for (String gene : allGenes) {
            Double expTValue = geneToExpTValue.get(gene);
            if (expTValue == null)
                continue;
            Set<String> fiPartners = fiToPartners.get(gene);
            fiPartners.retainAll(mutatedGenes);
            fu.printLine(gene + "\t" + fiPartners.size() + "\t" + expTValue);
        }
        fu.close();
    }
    
    /**
     * Analyzer the genes from differentially expressed genes, which are calculated from R,
     * and relationships with the clusters.
     * @throws Exception
     */
    @Test
    public void analyzeDiffExpGenesWithClusters() throws Exception {
        // Want to get FI partners for all mutated genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> fiToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        Set<String> allGenes  = new HashSet<String>(fiToPartners.keySet());
        System.out.println("Total genes in FI component: " + fiToPartners.size());
        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        System.out.println("Genes with two or more mutations: " + mutatedGenes.size());
        mutatedGenes.retainAll(fiToPartners.keySet());
        System.out.println("    in FI network component: " + mutatedGenes.size());
        Set<String> diffGenes = fu.loadInteractions(TCGA_GBM_DIR + "GBM_EXP_UP_GENES.txt");
        System.out.println("Hyper expression genes: " + diffGenes.size());
        diffGenes.retainAll(fiToPartners.keySet());
        System.out.println("    in FI network component: " + diffGenes.size());
        analyzeDiffGeneDist(allGenes, 
                            diffGenes,
                            mutatedGenes);
        // Check the distribution of these diff genes in the FI partners of mutated genes
        Set<String> mutatedGenePartners = new HashSet<String>();
        for (String gene : mutatedGenes) {
            Set<String> fiPartners = fiToPartners.get(gene);
            mutatedGenePartners.addAll(fiPartners);
        }
        mutatedGenePartners.removeAll(mutatedGenes);
        System.out.println("Total FI partners for mutated genes: " + mutatedGenePartners.size());
        analyzeDiffGeneDist(allGenes, 
                            diffGenes,
                            mutatedGenePartners);
        // Check the distribution of these diff genes in the core GBM cluster
        List<String> gbmCluster = getClusterForAnalysis(DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt",
                                                       "1.0");
        System.out.println("GBM cluster genes: " + gbmCluster.size());
        analyzeDiffGeneDist(allGenes, 
                            diffGenes,
                            gbmCluster);
        Set<String> gbmClusterPartners = new HashSet<String>();
        for (String gene : gbmCluster) {
            Set<String> fiPartners = fiToPartners.get(gene);
            gbmClusterPartners.addAll(fiPartners);
        }
        gbmClusterPartners.removeAll(gbmCluster);
        System.out.println("Total FI partners for gbm cluster: " + gbmClusterPartners.size());
        analyzeDiffGeneDist(allGenes, diffGenes, gbmClusterPartners);
    }
    
    private void analyzeDiffGeneDist(Collection<String> allGenes,
                                     Collection<String> diffGenes,
                                     Collection<String> targetGenes) throws Exception {
        Set<String> tmp = new HashSet<String>(diffGenes);
        tmp.retainAll(targetGenes);
        System.out.println("    shared genes: " + tmp.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(), 
                                                                    targetGenes.size(), 
                                                                    diffGenes.size(), 
                                                                    tmp.size());
        System.out.println("    pvalue: " + pvalue);
    }
    
    @Test
    public void checkExpGeneHierarchicalClusters() throws Exception {
        String fileName = TCGA_GBM_DIR + "HierarchicalClusterGeneExpUpDown.txt";
        HierarchicalCluster clustering = new HierarchicalCluster();
        List<HierarchicalClusterNode> clusters = clustering.loadHierarchicalClusters(fileName, 3.0);
        clustering.sortClustersBasedOnSizes(clusters);
        int index = 0;
        for (HierarchicalClusterNode cluster : clusters) {
            if (cluster.getIds().size() < 3)
                continue;
            System.out.println(index + "\t" + cluster.getIds().size());
            index ++;
        }
        // For gene expression
        Map<String, List<GeneExpressionValue>> sampleToExpValues = loadSampleToGeneExpValues();
        System.out.println("Total samples: " + sampleToExpValues.size());
        // Calculate average values for these clusters
        Map<String, double[]> sampleToClusterValues = new HashMap<String, double[]>();
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (String sample : sampleToExpValues.keySet()) {
            List<GeneExpressionValue> values = sampleToExpValues.get(sample);
            double[] clusterValues = new double[clusters.size()];
            index = 0;
            for (HierarchicalClusterNode cluster : clusters) {
                stat.clear();
                for (String gene : cluster.getIds()) {
                    for (GeneExpressionValue v : values) {
                        if (v.getGene().equals(gene)) {
                            stat.addValue(v.getZvalue());
                            break;
                        }
                    }
                }
                clusterValues[index] = stat.getMean();
                index ++;
            }
            sampleToClusterValues.put(sample, clusterValues);
        }
        // Output these values
        StringBuilder builder = new StringBuilder();
        for (String sample : sampleToClusterValues.keySet()) {
            double[] values = sampleToClusterValues.get(sample);
            builder.append(sample);
            for (double v : values)
                builder.append("\t").append(v);
            builder.append("\n");
        }
        System.out.println(builder.toString());
        // Calculate distance between samples
        Map<String, Double> sampleToDistance = new HashMap<String, Double>();
        List<String> sampleList = new ArrayList<String>(sampleToClusterValues.keySet());
        for (int i = 0; i < sampleList.size() - 1; i++) {
            String sample1 = sampleList.get(i);
            double[] values1 = sampleToClusterValues.get(sample1);
            for (int j = i + 1; j < sampleList.size(); j++) {
                String sample2 = sampleList.get(j);
                double[] values2 = sampleToClusterValues.get(sample2);
                // Want to get the distance
                double correlation = 1.0 - Utils.correlation(values1, 
                                                             values2,
                                                             values1.length);
                sampleToDistance.put(sample1 + "\t" + sample2, 
                                     correlation);
            }
        }
        // For a hierachical clustering
        System.out.println("hierarchical layouting...");
        List<HierarchicalClusterNode> sampleClusters = hierarchicalCluster(sampleList, 
                                                               sampleToDistance,
                                                               false);
        HierarchicalClusterNode root = sampleClusters.get(0);
        outputCluster(root);
    }
    
    /**
     * This method is used to analyze some topological relationships between highly expressed genes with 
     * mutated genes.
     * @throws Exception
     */
    @Test
    public void analyzeHyperExpGenesWithClusters() throws Exception {
        // Want to get FI partners for all mutated genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> fiToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        List<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        System.out.println("Genes with two or more mutations: " + mutatedGenes.size());
        mutatedGenes.retainAll(fiToPartners.keySet());
        System.out.println("    in FI network component: " + mutatedGenes.size());
        Set<String> hyperExpGenes = fu.loadInteractions(TCGA_GBM_DIR + "HighExpGenes5Percentile.txt");
        System.out.println("Hyper expression genes: " + hyperExpGenes.size());
        hyperExpGenes.retainAll(fiToPartners.keySet());
        System.out.println("    in FI network component: " + hyperExpGenes.size());
        Set<String> tmp = new HashSet<String>(hyperExpGenes);
        tmp.retainAll(mutatedGenes);
        System.out.println("Shared genes: " + tmp.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(fiToPartners.size(), 
                                                                    mutatedGenes.size(), 
                                                                    hyperExpGenes.size(), 
                                                                    tmp.size());
        System.out.println("    pvalue: " + pvalue);
        Set<String> mutatedGenePartners = new HashSet<String>();
        for (String gene : mutatedGenes) {
            Set<String> fiPartners = fiToPartners.get(gene);
            mutatedGenePartners.addAll(fiPartners);
        }
        mutatedGenePartners.removeAll(mutatedGenes);
        System.out.println("Total FI partners for mutated genes: " + mutatedGenePartners.size());
        tmp = new HashSet<String>(hyperExpGenes);
        tmp.retainAll(mutatedGenePartners);
        System.out.println("Hyper genes in FI partners: " + tmp.size());
        pvalue = MathUtilities.calculateHypergeometricPValue(fiToPartners.size(), 
                                                             mutatedGenePartners.size(), 
                                                             hyperExpGenes.size(), 
                                                             tmp.size());
        System.out.println("    pvalue: " + pvalue);
    }
    
    private List<String> getControlSamples() {
        String[] names = new String[] {
                "TCGA-06-0673-11",
                "TCGA-06-0675-11",
                "TCGA-06-0676-11",
                "TCGA-06-0678-11",
                "TCGA-06-0680-11",
                "TCGA-06-0681-11",
                "TCGA-08-0623-11",
                "TCGA-08-0625-11",
                "TCGA-08-0626-11",
                "TCGA-08-0627-11",
                "TCGA-07-0249-20"
        };
        return Arrays.asList(names);
    }
    
    /**
     * Generate a type for R multtest.
     * @throws IOException
     */
    @Test
    public void generateTypes() throws IOException {
        String fileName = TCGA_GBM_DIR + "HG-U133A_Gene_Exp_table.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        fu.close();
        List<String> controls = getControlSamples();
        int total = 0;
        int controlNumber = 0;
        for (int i = 1; i < headers.length; i++) {
            String header = headers[i];
            boolean isControl = isControlSample(header, controls);
            if (isControl)
                controlNumber ++;
            System.out.print((isControl ? "0" : "1") + " ");
            total ++;
        }
        System.out.println();
        System.out.println("Total: " + total);
        System.out.println("Control: " + controlNumber);
    }
    
    private boolean isControlSample(String sample, 
                                    List<String> controls) {
        for (String control : controls) {
            if (sample.startsWith(control)) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * This method is used to choose the first one percentile genes across all samples.
     * @throws IOException
     */
    @Test
    public void chooseHyperExpGenes() throws IOException {
        String fileName = TCGA_GBM_DIR + "HG-U133A_Gene_Exp_table.txt";
        // Load all genes
        final Map<String, List<Double>> geneToValues = new HashMap<String, List<Double>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            List<Double> values = new ArrayList<Double>();
            for (int i = 1; i < tokens.length; i++) {
                Double value = new Double(tokens[i]);
                values.add(value);
            }
            geneToValues.put(gene, values);
        }
        fu.close();
        Set<String> sharedGenes = new HashSet<String>();
        int firstPercentile = (int) (geneToValues.size() * 0.05);
        List<String> controls = getControlSamples();
        Set<String> controlGenes = new HashSet<String>();
        for (int i = 1; i < headers.length; i++) {
            String sampleName = headers[i];
            // Check if this sample is a control
            boolean isControl = isControlSample(headers[i], 
                                                controls);
            Set<String> workingGenes = null;
            if (isControl)
                workingGenes = controlGenes;
            else
                workingGenes = sharedGenes;
            List<String> geneList = new ArrayList<String>(geneToValues.keySet());
            final int index = i;
            Collections.sort(geneList, new Comparator<String>() {
                public int compare(String gene1, String gene2) {
                    List<Double> values1 = geneToValues.get(gene1);
                    List<Double> values2 = geneToValues.get(gene2);
                    Double value1 = values1.get(index - 1);
                    Double value2 = values2.get(index - 1);
                    return value2.compareTo(value1);
                }
            });
            // Pick up the first one-percentile genes
            List<String> tmp = new ArrayList<String>();
            for (int j = 0; j < firstPercentile; j++) {
                tmp.add(geneList.get(j));
            }
            if (workingGenes.size() == 0)
                workingGenes.addAll(tmp);
            else
                workingGenes.retainAll(tmp);
        }
        System.out.println("Total shared genes: " + sharedGenes.size());
        System.out.println("Total control genes: " + controlGenes.size());
        sharedGenes.removeAll(controlGenes); // More like house-keep genes
        System.out.println("Removing genes from control: " + sharedGenes.size());
        for (String gene : sharedGenes)
            System.out.println(gene);
        String outFileName = TCGA_GBM_DIR + "HighExpGenes5Percentile.txt";
        fu.setOutput(outFileName);
        List<String> list = new ArrayList<String>(sharedGenes);
        Collections.sort(list);
        for (String gene : list)
            fu.printLine(gene);
        fu.close();
    }
    
    /**
     * Convert a list of expression value into a matrix
     * @throws IOException
     */
    @Test
    public void transformGeneExpDataFile() throws IOException {
        String dirName = "datasets/TCGA/GBM/";
        String inFileName = dirName + "HG-U133A_Gene_Exp.txt";
        String outFileName = dirName + "HG-U133A_Gene_Exp_table.txt";
        fu.setInput(inFileName);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        String line = fu.readLine();
        // Get all samples
        List<String> samples = new ArrayList<String>();
        Map<String, List<String>> geneToValues = new HashMap<String, List<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!samples.contains(tokens[0]))
                samples.add(tokens[0]);
            List<String> list = geneToValues.get(tokens[1]);
            if (list == null) {
                list = new ArrayList<String>();
                geneToValues.put(tokens[1], list);
            }
            list.add(tokens[2]);
        }
        // Want to print out 
        StringBuilder builder = new StringBuilder();
        // for headers
        builder.append("Gene");
        for (String sample : samples) {
            builder.append("\t").append(sample);
        }
        outFu.printLine(builder.toString());
        builder.setLength(0);
        List<String> geneList = new ArrayList<String>(geneToValues.keySet());
        Collections.sort(geneList);
        for (String gene : geneList) {
            List<String> values = geneToValues.get(gene);
            if (values.size() < samples.size()) {
                System.err.println(gene + " has not enough data points!");
                continue;
            }
            builder.append(gene);
            for (String value : values) {
                builder.append("\t").append(value);
            }
            outFu.printLine(builder.toString());
            builder.setLength(0);
        }
        outFu.close();
        fu.close();
    }
    
    /**
     * This method is used to check the graph component in the mutated genes.
     * @throws Exception
     */
    @Test
    public void checkGraphComponentsInMutatedGenes() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiIds = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Integer> geneToMutation = loadNatureGBMMutationTable();
        geneToMutation.keySet().retainAll(fiIds);
        filterGenesBasedOnMutationNumber(geneToMutation, 2);
        Set<String> mutatedGenes = new HashSet<String>(geneToMutation.keySet());
        //Set<String> mutatedGenes = loadAllNatureGBM601Genes();
        int totalMutatedGenes = mutatedGenes.size();
        System.out.println("Total mutated genes: " + totalMutatedGenes);
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        Set<String> fiCopy = new HashSet<String>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            if (mutatedGenes.contains(gene1) &&
                mutatedGenes.contains(gene2))
                fiCopy.add(fi);
        }
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fiCopy);
        SimpleGraph<String, DefaultEdge> graph = (SimpleGraph<String, DefaultEdge>) JGraphTUtilities.createGraph(ids, fiCopy);
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<String, DefaultEdge>(graph);
        List<Set<String>> components = inspector.connectedSets();
        System.out.println("Total components: " + components.size());
        // Search the biggest component
        Set<String> biggest = new HashSet<String>();
        for (Set<String> tmp : components) {
            if (tmp.size() > biggest.size())
                biggest = tmp;
        }
        System.out.println("Biggest component: " + biggest.size() + " (" + (double)biggest.size() / mutatedGenes.size() + ")");
    }
    
    @Test
    public void checkPreSelectedGenesInFINetwork() throws Exception {
        Set<String> allGenes = loadAllNatureGBM601Genes();
        System.out.println("Total GBM genes: " + allGenes.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        allGenes.retainAll(fiGenes);
        System.out.println("Genes in FI network: " + allGenes.size());
    }
    
    /**
     * This method is used to check FI partners for all 601 resquenced genes.
     * @throws Exception
     */
    @Test
    public void checkFIPartnersForAllGenes() throws Exception {
        Set<String> allGenes = loadAllNatureGBM601Genes();
        System.out.println("Total GBM genes: " + allGenes.size());
        // Check connection degrees
        // Want to check the biggest components only
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToFIPartners = bfs.generateIdToPartnersMap(fis);
        allGenes.retainAll(geneToFIPartners.keySet());
        int totalGBMGenes = allGenes.size();
        System.out.println("Total GBM genes in FI network: " + totalGBMGenes);
        List<String> allGeneList = new ArrayList<String>(allGenes);
        Collections.sort(allGeneList);
        int total = geneToFIPartners.size();
        System.out.println("Total genes: " + total);
        System.out.println("Gene\tFIPartners\tGBMFIPartners\tP-value");
        for (String gene : allGeneList) {
            Set<String> partners = geneToFIPartners.get(gene);
            int fiPartners = partners.size();
            partners.retainAll(allGenes);
            int GBMFIPartners = partners.size();
            double pvalue = MathUtilities.calculateHypergeometricPValue(total - 1, // Remove itself
                                                                        fiPartners, 
                                                                        totalGBMGenes - 1, 
                                                                        GBMFIPartners);
            System.out.println(gene + "\t" + fiPartners + "\t" +
                               GBMFIPartners + "\t" + pvalue);
        }
    }
    
    /**
     * This method is used to check FI partners for mutated genes from NatureGBM.
     * @throws IOException
     */
    @Test
    public void checkFIPartnersForMutatedGenes() throws Exception {
        Set<String> allGBMGenes = loadAllNatureGBM601Genes();
        Map<String, Integer> geneToMutation = loadNatureGBMMutationTable();
        //filterGenesBasedOnMutationNumber(geneToMutation, 2);
        int totalMutated = geneToMutation.size();
        System.out.println("Total mutated genes: " + totalMutated);
        Set<String> mutatedGeneSet = new HashSet<String>(geneToMutation.keySet());
        // Check connection degrees
        // Want to check the biggest components only
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToFIPartners = bfs.generateIdToPartnersMap(fis);
        // Just want to see genes in the FI network
        mutatedGeneSet.retainAll(geneToFIPartners.keySet()); 
        allGBMGenes.retainAll(geneToFIPartners.keySet());
        totalMutated = mutatedGeneSet.size();
        int totalGBMGenes = allGBMGenes.size();
        System.out.println("Total mutated genes in FI network: " + totalMutated);
        System.out.println("Total GBM genes in FI network: " + totalGBMGenes);
        List<String> mutatedGenesList = new ArrayList<String>(mutatedGeneSet);
        Collections.sort(mutatedGenesList);
        int total = geneToFIPartners.size();
        System.out.println("Total genes: " + total);
        System.out.println("Gene\tFIPartners\tMutatedFIPartners\tP-value");
        for (String gene : mutatedGenesList) {
            Integer number = geneToMutation.get(gene);
            if (number < 3)
                continue;
            Set<String> partners = geneToFIPartners.get(gene);
            partners.retainAll(allGBMGenes); // Want to calculate in the GBM gene set
            int fiPartners = partners.size();
            partners.retainAll(mutatedGeneSet);
            int mutatedFIPartners = partners.size();
            double pvalue = MathUtilities.calculateHypergeometricPValue(totalGBMGenes - 1, // Remove itself
                                                                        fiPartners, 
                                                                        totalMutated - 1, 
                                                                        mutatedFIPartners);
            System.out.println(gene + "\t" + fiPartners + "\t" +
                               mutatedFIPartners + "\t" + pvalue);
        }
    }
    
    private void filterGenesBasedOnMutationNumber(Map<String, Integer> geneToMutation, int cutoff) {
        for (Iterator<String> it = geneToMutation.keySet().iterator(); it.hasNext();) {
            String gene = it.next();
            Integer number = geneToMutation.get(gene);
            if (number < cutoff)
                it.remove();
        }
    }
    
    @Test
    public void checkMutationNumbers() throws IOException {
        final Map<String, Integer> geneToMutation = loadNatureGBMMutationTable();
        List<String> genes = new ArrayList<String>(geneToMutation.keySet());
        Collections.sort(genes, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Integer number1 = geneToMutation.get(gene1);
                Integer number2 = geneToMutation.get(gene2);
                return number2.compareTo(number1);
            }
        });
        for (String gene : genes) {
            Integer number = geneToMutation.get(gene);
            System.out.println(gene + "=" + number);
        }
    }
    
    /**
     * This method is used to generate a pair-wise relationships between degrees and
     * mutation numbers.
     * @throws IOException
     */
    @Test
    public void generateMutationNumbersVsDegrees() throws Exception {
        List<String> samples = loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = getSampleToAlteredGenes(samples);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToAlteredGenes);
        String[] intFileNames = new String[] {
                R3Constants.RESULT_DIR + "FIsInGene_041709.txt", // All FIs
                R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt", // Pathway FIs
                R3Constants.RESULT_DIR + "FIsInGene_Predicted_041709.txt" // Predicted FIs
        };
        String[] outFileNames = new String[] {
                TCGA_GBM_DIR + "NatureGBNAlterationToDegree_All.txt",
                TCGA_GBM_DIR + "NatureGBNAlterationToDegree_Pathway.txt",
                TCGA_GBM_DIR + "NatureGBNAlterationToDegree_Predicted.txt"
        };
        for (int i = 0; i < intFileNames.length; i++) {
            String intFileName = intFileNames[i];
            String outFileName = outFileNames[i];
            Set<String> fis = fu.loadInteractions(intFileName);
            BreadthFirstSearch bfs = new BreadthFirstSearch();
            //bfs.setInteractions(interactions);
            Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
            fu.setOutput(outFileName);
            for (String gene : geneToSamples.keySet()) {
                Set<String> samples1 = geneToSamples.get(gene);
                Integer mutation = samples1.size();
                Set<String> partners = geneToPartners.get(gene);
                if (partners == null)
                    continue;
                String line = gene + "\t" + mutation + "\t" + partners.size();
                System.out.println(line);
                fu.printLine(line);
            }
            fu.close();
        }
    }
    
}
