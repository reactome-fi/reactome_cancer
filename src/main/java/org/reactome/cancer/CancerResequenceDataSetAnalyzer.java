/*
 * Created on Sep 10, 2008
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.commons.math.stat.inference.TestUtils;
import org.hibernate.Session;
import org.jgrapht.Graph;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.hibernate.FISourceTypeReader;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.r3.UCSCDataAnalyzer;
import org.reactome.r3.cluster.DistanceCalculator;
import org.reactome.r3.cluster.HierarchicalCluster;
import org.reactome.r3.cluster.HierarchicalClusterNode;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.BreadthFirstSearch.Edge;
import org.reactome.r3.graph.BreadthFirstSearch.TreeNode;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.NetworkBuilderForGeneSet;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to analyze the data set related to Pancreatic cancer.
 * @author wgm
 *
 */
public class CancerResequenceDataSetAnalyzer implements CancerConstants {
    //private final int HUB_CUTOFF = 127; // Proteins with degree greater than 127 are hubs (exclusive)
    private final int HUB_CUTOFF = 219;
    protected final String DIR_NAME = "datasets/PancreaticCancers/";
    private final String GENE_FI_FILE_NAME = R3Constants.GENE_FI_FILE_NAME;
    private final String GENE_FI_BIG_COMP_FILE_NAME = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
    // Check with names in upper case
    //private final String GENE_FI_FILE_NAME = R3Constants.RESULT_DIR + "FI73InGeneUpperCase_111208.txt";
    //private final String GENE_FI_BIG_COMP_FILE_NAME = R3Constants.RESULT_DIR + "FI73InGeneUpperCase_111208_BigComp.txt";
    protected FileUtility fu = new FileUtility();
    public static final String OVARIAN_DIR_NAME = R3Constants.OVARIAN_DIR_NAME;
    protected final String TCGA_GBM_DIR_NAME = R3Constants.GBM_DIR;
    protected final String OVARIAN_RESULT_DIR = OVARIAN_DIR_NAME + "results/";
    
    public static void main(String[] args) throws Exception {
        CancerResequenceDataSetAnalyzer analyzer = new CancerResequenceDataSetAnalyzer();
        analyzer.permutationTestForMinSpanBasedOnCluster();
    }
    
    public CancerResequenceDataSetAnalyzer() {
    }    
    
    /**
     * Do some quick stat analysis
     * @param sampleToGenes
     * @throws Exception
     */
    public void summarizeSampleToGenes(Map<String, Set<String>> sampleToGenes)  {
        SummaryStatistics summary = new SummaryStatistics();
        for (Set<String> genes : sampleToGenes.values())
            summary.addValue(genes.size());
        System.out.println("Mean: " + summary.getMean());
        System.out.println("Max: " + summary.getMax());
        System.out.println("Min: " + summary.getMin());
        System.out.println("SD: " + summary.getStandardDeviation());
    }
    
    /**
     * This method is used to convert a FI file in name to another FI file in UniProt.
     * @throws Exception
     */
    @Test
    public void convertFIInNameToFIInUniProt() throws Exception {
        String[] inFiles = new String[] {
                "FIsInCluster17ForScienceGBMViaHierarchical102908.txt",
                "FIsInCluster57ForNatureGBMViaHierarchical102908.txt",
                "FIsInCluster30ForSciencePancreaticViaHierarchical102908.txt"
        };
        String[] outFiles = new String[] {
                "FIsInCluster17ForScienceGBMViaHierarchical102908InUniProt.txt",
                "FIsInCluster57ForNatureGBMViaHierarchical102908InUniProt.txt",
                "FIsInCluster30ForSciencePancreaticViaHierarchical102908InUniProt.txt"
        };
        HibernateFIReader fiReader = new HibernateFIReader();
        Map<String, Set<String>> nameToAccesses = fiReader.generateProteinNameToAccession();
        FileUtility outFu = new FileUtility();
        for (int i = 0; i < inFiles.length; i++) {
            String inFile = inFiles[i];
            Set<String> fis = fu.loadInteractions(DIR_NAME + inFile);
            String outFile = outFiles[i];
            outFu.setOutput(DIR_NAME + outFile);
            for (String fi : fis) {
                int index = fi.indexOf("\t");
                String name1 = fi.substring(0, index);
                String name2 = fi.substring(index + 1);
                Set<String> accessions1 = nameToAccesses.get(name1);
                accessions1 = cleanUpAccessions(accessions1);
                Set<String> accessions2 = nameToAccesses.get(name2);
                accessions2 = cleanUpAccessions(accessions2);
                // Permutate to create FIs
                for (String acc1 : accessions1) {
                    for (String acc2 : accessions2) {
                        if (acc1.equals(acc2))
                            continue; // Don't want to have self-interaction
                        outFu.printLine(acc1 + "\t" + acc2);
                    }
                }
            }
            outFu.close();
        }
    }
    
    /**
     * A helper method to remove non-uniprot ids and splice forms.
     * @param accessions
     */
    private Set<String> cleanUpAccessions(Set<String> accessions) {
        Set<String> rtn = new HashSet<String>();
        for (Iterator<String> it = accessions.iterator(); it.hasNext();) {
            String accession = it.next();
            if (accession.matches("(^\\d)(\\w)*")) {
                //System.out.println(accession);
                //it.remove();
                continue;
            }
            int index = accession.indexOf("-");
            if (index > 0)
                rtn.add(accession.substring(0, index));
            else
                rtn.add(accession);
        }
        return rtn;
    }
    
    @Test
    public void checkFIsInGenes() throws Exception {
        String geneList = "PSG1, VNN1, TNFRSF10C, NOVA2, ZNF800, " +
        "NIPBL, BBC3, HNRNPUL2, NUCB2, LIPN, ERGIC2, MLLT4, " +
        "FBXL20, CHAD, LOC641367, SLC4A8, FEN1, ZNF426, TARS, C12orf44";
        String[] tmp = geneList.split(", ");
        List<String> names = Arrays.asList(tmp);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fisInGenes = InteractionUtilities.getFIs(names, fis);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> shared = InteractionUtilities.getShared(allGenes, names);
        System.out.println("Genes in the FI network: " + shared.size());
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        System.out.println("Total genes for query: " + names.size());
        System.out.println("Total FIs for genes: " + fisInGenes.size());
        System.out.println("Total genes having FIs: " + genesInFIs.size());
    }
    
    /**
     * This method is used to check cluster and cutoff relationships.
     * @throws IOException
     */
    @Test
    public void checkClusterAndCutoff() throws IOException {
        String fileName = DIR_NAME + "HierarchicalClusters102908.txt";
        fu.setInput(fileName);
        Map<Double, Integer> scoreToNumber = new HashMap<Double, Integer>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("datasets")) {
                // A new file
                System.out.println(line);
            }
            else if (line.length() == 0) {
                outputClusterScoreToSize(scoreToNumber);
                scoreToNumber.clear();
            }
            else {
                int index = line.indexOf(":");
                if (index == -1)
                    continue; // A single node
                Double score = Double.parseDouble(line.substring(0, index));
                index = line.lastIndexOf("(");
                Integer number = Integer.parseInt(line.substring(index + 1, 
                                                                 line.length() - 1));
                Integer existed = scoreToNumber.get(score);
                if (existed == null || existed < number) {
                    scoreToNumber.put(score, number);
                }
            }
        }
        fu.close();
        outputClusterScoreToSize(scoreToNumber);
    }
    
    private void outputClusterScoreToSize(Map<Double, Integer> scoreToNumber) {
        List<Double> scoreList = new ArrayList<Double>(scoreToNumber.keySet());
        Collections.sort(scoreList);
        int preSize = 0;
        System.out.println("score\tsize");
        for (Double score : scoreList) {
            Integer number = scoreToNumber.get(score);
            if (number > preSize)
                preSize = number;
            System.out.println(score + "\t" + preSize);
        }
    }
    
    /**
     * This method is used to check the linked genes for subnetwork.
     * @throws IOException
     */
    @Test
    public void checkLinkerGenes() throws IOException {
        Set<String> allNatureMutatedGenes = loadNatureGlioblastomaGenes();
        Set<String> allScienceMutatedGenes = loadAllScienceGBMMutatedGenes();
        Set<String> allMutatedGenes = new HashSet<String>();
        allMutatedGenes.addAll(allNatureMutatedGenes);
        allMutatedGenes.addAll(allScienceMutatedGenes);
        
        String[] fiFileNames = new String[] {
                "FIsInCluster17ForScienceGBMViaHierarchical102908.txt",
                "FIsInCluster30ForSciencePancreaticViaHierarchical102908.txt",
                "FIsInCluster57ForNatureGBMViaHierarchical102908.txt"
        };
        String[] clusterFileNames = new String[] {
                "ClusterListForScienceGBM102908.txt",
                "ClusterListForSciencePancreatic102908.txt",
                "ClusterListForNatureGBM102908.txt"
        };
        String[] cutoffs = new String[] {
                "2.94",
                "3.15",
                "3.18"
        };
        for (int i = 0; i < fiFileNames.length; i++) {
            System.out.println(fiFileNames[i]);
            // Want to check how many genes have mutation (this should be only one mutation from
            // the whole screen since we have picked up genes having at least two mutations)
            String fiFileName = DIR_NAME + fiFileNames[i];
            String clusterFileName = DIR_NAME + clusterFileNames[i];
            String cutoff = cutoffs[i];
            
//          String fiFileName = DIR_NAME + "FIsInCluster52ForNatureGBMViaHierarchical.txt";
//          String clusterFileName = DIR_NAME + "ClusterListForGlioblastoma.txt";
//          String cutoff = "3.06";
            List<String> clusterGenes = getClusterForAnalysis(clusterFileName,
                                                              cutoff);
            System.out.println("Total cluster genes: " + clusterGenes.size());
            Set<String> fis = fu.loadInteractions(fiFileName);
            Set<String> allFiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
            // Find how many genes are linker genes
            System.out.println("Total FI genes: " + allFiGenes.size());
            allFiGenes.removeAll(clusterGenes);
            System.out.println("Total linker genes: " + allFiGenes.size());
            System.out.println("Total all mutated genes: " + allMutatedGenes.size());
            allFiGenes.retainAll(allMutatedGenes);
            System.out.println("Total linker genes having mutation: " + allFiGenes);
        }
    }
    
    protected List<String> getClusterForAnalysis(String fileName, 
                                                String cutoff) throws IOException {
        Map<String, List<String>> clusters = loadClusters(fileName);
        // Get the cluster
        List<String> cluster = null;
        for (String key : clusters.keySet()) {
            if (key.contains(cutoff)) {
                cluster = clusters.get(key);
                break;
            }
        }
        return cluster;
    }

    
    /**
     * This method is used to analyze overlapping GBM genes published in Nature and Science.
     * @throws IOException
     */
    @Test
    public void checkSharedGBMGenesFromNatureAndScience() throws IOException {
        List<String> scienceGenes = loadScienceGBMGenes();
        List<String> natureGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        System.out.println("Total genes from Science: " + scienceGenes.size());
        System.out.println("Total genes from Nature: " + natureGenes.size());
        natureGenes.retainAll(scienceGenes);
        System.out.println("Shared genes: " + natureGenes.size());
        System.out.println("Shared genes: " + natureGenes);
        // Check shared genes in the clusters
        List<String> scienceCluster = getClusterForAnalysis(DIR_NAME + "ClusterListForScienceGBM042109.txt", 
                                                            "3.0625");
        List<String> natureCluster = getClusterForAnalysis(DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt",
                                                           "1.0");
        // Comparing clusters from average path
        System.out.println("Compare clusters between Nature and Science GBMs:");
        compareTwoGeneLists(natureCluster,
                            scienceCluster);
    }
    
    private List<String> getSharedGBMGenesFromNatureAndScience() throws IOException {
        List<String> scienceGenes = loadScienceGBMGenes();
        List<String> natureGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        natureGenes.retainAll(scienceGenes);
        return natureGenes;
    }
    
    /**
     * This method is used to check degreee for cancer genes in clusters.
     * @throws IOException
     */
    @Test
    public void checkProteinDegrees() throws IOException {
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        List<String> cancerGenes = getSharedGBMGenesFromNatureAndScience();
        for (String gene : cancerGenes) {
            Integer degree = proteinToDegree.get(gene);
            System.out.println(gene + "\t" + degree);
        }
    }
    
    private void compareTwoGeneLists(List<String> list1,
                                     List<String> list2) {
        System.out.println("First list: " + list1.size());
        System.out.println("Second list: " + list2.size());
        List<String> copy1 = new ArrayList<String>(list1);
        copy1.retainAll(list2);
        System.out.println("Shared genes: " + copy1.size());
    }
    
    @Test
    public void permutationTestForMinSpanBasedOnClusterBasedOnDegree() throws IOException {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        // These statements are used for degree based permutation
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        // Want to permutated pre-selected GBM genes
//        Set<String> allGBMGenes = loadAllNatureGBM601Genes();
//        proteinToDegree.keySet().retainAll(allGBMGenes);
        System.out.println("Total sampling genes: " + proteinToDegree.size());
 
//        List<Set<String>> bins = sortAllGenesBasedOnDegrees(proteinToDegree);
        //        Set<String> hubProteins = new HashSet<String>();
//        Set<String> nonHubProteins = new HashSet<String>();
//        sortProteinsForHubs(proteinToDegree, hubProteins, nonHubProteins);
//        // For hierarchical clustering
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        int permutationNumber = 1000;
        double[] cutoffs = new double[] {
//                3.18,
//                2.94,
//                3.15
                //1.0
                2.9231
        };
        String clusterFileNames[] = new String[] {
//                DIR_NAME + "ClusterListForNatureGBM102908.txt",
//                DIR_NAME + "ClusterListForScienceGBM102908.txt",
//                DIR_NAME + "ClusterListForSciencePancreatic102908.txt"
//                DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt"
                OVARIAN_RESULT_DIR + "ClusterList051309.txt"
        };
        // Permutation test based on bins generated in method sortAllGenesBasedOnDegrees(Map<String, Ineger>)
//        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeBin_NatureGBM_1000_111408.txt",
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeBin_ScienceGBM_1000_111408.txt",
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeBin_SciencePancreatic_1000_111408.txt"
//        };
        // Permutation test based on bins generated in method sortAllGenesBasedOnDegreesInCancerGenes
        // (Map<String, Integer>, List<String>)
        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeSelfBin_NatureGBM_1000_111508.txt",
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeSelfBin_ScienceGBM_1000_111508.txt",
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeSelfBin_SciencePancreatic_1000_111508.txt"
//                DIR_NAME + "PermutationTestForMiniSpanViaClusterBasedOnDegreeSelfBinPre_NatureGBM_Cut_1_1000_111508.txt"
                OVARIAN_RESULT_DIR + "PermutationTestForMiniSpanViaClusterBasedOnDegreeSelfBin_1000_051309.txt",
        };        
        HierarchicalCluster hclust = new HierarchicalCluster();
        for (int i = 0; i < cutoffs.length; i++) {
            // Need to get genes
            String clusterFileName = clusterFileNames[i];
            double cutoff = cutoffs[i];
            List<String> genes = getClusterForAnalysis(clusterFileName, cutoff + "");
            List<Set<String>> bins = sortAllGenesBasedOnDegreesInCancerGenes(proteinToDegree, 
                                                                             genes);
            List<Integer> sizeInBins = calculateNumbersInBins(genes, 
                                                              bins);
//            // Find how many in hubs, how many in non-hubs
//            List<String> hubInGenes = new ArrayList<String>(genes);
//            hubInGenes.retainAll(hubProteins);
//            List<String> nonHubInGenes = new ArrayList<String>(genes);
//            nonHubInGenes.removeAll(hubInGenes);
            System.out.println(clusterFileName);
            System.out.println("Starting permutation...");
            List<Integer> weights = new ArrayList<Integer>();
            for (int j = 0; j < permutationNumber; j ++) {
                //long time1 = System.currentTimeMillis();
//                Set<String> sampled = new HashSet<String>();
//                Set<String> sampled1 = randomSampling(hubProteins, hubInGenes.size());
//                Set<String> sampled2 = randomSampling(nonHubProteins, nonHubInGenes.size());
//                sampled.addAll(sampled1);
//                sampled.addAll(sampled2);
                Set<String> sampled = samplingBasedOnBins(sizeInBins, bins);
                List<String> geneList = new ArrayList<String>(sampled);
                Map<String, Integer> pairToDistance = calculateDistances(geneList);
                List<HierarchicalClusterNode> newClusters = hierarchicalCluster(geneList, 
                                                                    pairToDistance, 
                                                                    false);
                int totalWeight = 0;
                // Calculate total cluster weight
                HierarchicalClusterNode newFirstNode = newClusters.get(0);
                Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
                hclust.grepAllClusters(newFirstNode, allNodes);
                //for (ClusterNode node : allNodes)
                //    totalWeight += node.pathDistance;
                //System.out.println("Total weight: " + totalWeight);
                // These are used to calculate the actual path
                Map<String, List<String>> pairToPath = bfs.generateShortestPath(geneList, 
                                                                                nodeToEdges);
                Set<String> spanningFIs = generateSpanFromClusters(allNodes,
                                                                   pairToPath);
                weights.add(spanningFIs.size());
            }
            Collections.sort(weights);
            String outFileName = outFileNames[i];
            fu.setOutput(outFileName);
            for (Integer ratio : weights)
                fu.printLine(ratio + "");
            fu.close();
        }
    }
    
    /**
     * Permutation test for pre-defined list for total weights of MSTs.
     * @throws IOException
     */
    @Test
    public void permutationTestForMinSpanBasedOnClusterPreList() throws IOException {
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        // Prepare proteins
        Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(fis);
        // Load genes
        Set<String> preList = loadAllNatureGBM601Genes();
        preList.retainAll(proteins);
        
        int permutationNumber = 1000;
        int[] sampleSizes = new int[] {
                57,
                //17,
                //30
                //21
        };
        String[] fileNames = new String[] {
                DIR_NAME + "PermutationTestForMiniSpanViaClusterPreList_57_1000_112108_1.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_17_1000_103008.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_30_1000_103008.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_21_1000_110308.txt",
        };
        HierarchicalCluster hclust = new HierarchicalCluster();
        for (int j = 0; j < sampleSizes.length; j++) {
            int sampleSize = sampleSizes[j];
            String fileName = fileNames[j];
            List<Integer> weights = new ArrayList<Integer>();
            for (int i = 0; i < permutationNumber; i++) {
                //long time1 = System.currentTimeMillis();
                Set<String> sample = randomSampling(preList, sampleSize);
                List<String> geneList = new ArrayList<String>(sample);
                Map<String, Integer> pairToDistance = calculateDistances(geneList);
                List<HierarchicalClusterNode> newClusters = hierarchicalCluster(geneList, 
                                                                    pairToDistance, 
                                                                    true);
                // Calculate total cluster weight
                HierarchicalClusterNode newFirstNode = newClusters.get(0);
                Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
                hclust.grepAllClusters(newFirstNode, allNodes);
                //for (ClusterNode node : allNodes)
                //    totalWeight += node.pathDistance;
                //System.out.println("Total weight: " + totalWeight);
                // These are used to calculate the actual path
                Map<String, List<String>> pairToPath = bfs.generateShortestPath(geneList, 
                                                                                nodeToEdges);
                Set<String> spanningFIs = generateSpanFromClusters(allNodes,
                                                                   pairToPath);
                weights.add(spanningFIs.size());
                //System.out.println(i + ": " + "Total FIs: " + spanningFIs.size());
                //long time2 = System.currentTimeMillis();
                //System.out.println("Time: " + (time2 - time1));
            }
            Collections.sort(weights);
            fu.setOutput(fileName);
            for (Integer weight : weights)
                fu.printLine(weight + "");
            fu.close();
        }
    }
    
    /**
     * This method is used to do permutation test for mini-span sub graph for a list of
     * genes randomly generated from the network.
     * @throws IOException
     */
    @Test
    public void permutationTestForMinSpanBasedOnCluster() throws IOException {
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        // Prepare proteins
        Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> proteinList = new ArrayList<String>(proteins);
        int permutationNumber = 1000;
        
        int[] sampleSizes = new int[] {
                //57,
                //17,
                //30
                //21
                11
        };
        String[] fileNames = new String[] {
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_57_1000_103008.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_17_1000_103008.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_30_1000_103008.txt",
                //DIR_NAME + "PermutationTestForMiniSpanViaCluster_21_1000_110308.txt",
                DIR_NAME + "PermutationTestForMiniSpanViaCluster_11_1000_041709.txt",
        };
        HierarchicalCluster clusteringHelper = new HierarchicalCluster();
        for (int j = 0; j < sampleSizes.length; j++) {
            int sampleSize = sampleSizes[j];
            String fileName = fileNames[j];
            List<Integer> weights = new ArrayList<Integer>();
            for (int i = 0; i < permutationNumber; i++) {
                //long time1 = System.currentTimeMillis();
                Set<String> sample = randomSampling(proteinList, sampleSize);
                List<String> geneList = new ArrayList<String>(sample);
                Map<String, Integer> pairToDistance = calculateDistances(geneList);
                List<HierarchicalClusterNode> newClusters = hierarchicalCluster(geneList, 
                                                                    pairToDistance, 
                                                                    true);
                int totalWeight = 0;
                // Calculate total cluster weight
                HierarchicalClusterNode newFirstNode = newClusters.get(0);
                Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
                clusteringHelper.grepAllClusters(newFirstNode, allNodes);
                //for (ClusterNode node : allNodes)
                //    totalWeight += node.pathDistance;
                //System.out.println("Total weight: " + totalWeight);
                // These are used to calculate the actual path
                Map<String, List<String>> pairToPath = bfs.generateShortestPath(geneList, 
                                                                                nodeToEdges);
                Set<String> spanningFIs = generateSpanFromClusters(allNodes,
                                                                   pairToPath);
                weights.add(spanningFIs.size());
                //System.out.println(i + ": " + "Total FIs: " + spanningFIs.size());
                //long time2 = System.currentTimeMillis();
                //System.out.println("Time: " + (time2 - time1));
            }
            Collections.sort(weights);
            fu.setOutput(fileName);
            for (Integer weight : weights)
                fu.printLine(weight + "");
            fu.close();
        }
    }
    
    /**
     * This method is used to calculcate the path for the mini-span network.
     * @param geneNames
     * @return
     * @throws Exception
     */
    public int calculatePathForMinispanNetwork(List<String> geneList) throws Exception {
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        Map<String, Integer> pairToDistance = calculateDistances(geneList);
        List<HierarchicalClusterNode> newClusters = hierarchicalCluster(geneList, 
                                                            pairToDistance, 
                                                            true);
        int totalWeight = 0;
        // Calculate total cluster weight
        HierarchicalClusterNode newFirstNode = newClusters.get(0);
        HierarchicalCluster cluster = new HierarchicalCluster();
        Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
        cluster.grepAllClusters(newFirstNode, allNodes);
        //for (ClusterNode node : allNodes)
        //    totalWeight += node.pathDistance;
        //System.out.println("Total weight: " + totalWeight);
        // These are used to calculate the actual path
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(geneList, 
                                                                        nodeToEdges);
        Set<String> spanningFIs = generateSpanFromClusters(allNodes,
                                                           pairToPath);
        return spanningFIs.size();
    }
    
    @Test
    public void checkTwoDistancesMethods() throws IOException {
        String clusterFileName = DIR_NAME + "ClusterListForNatureGBM102908.txt";
        String cutoff = "3.18";
        List<String> clusterGenes = getClusterForAnalysis(clusterFileName, cutoff);
        System.out.println("Total genes: " + clusterGenes.size());
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        Map<String, Integer> pairToDistance1 = calculateDistances(clusterGenes,
                                                                  bfs, 
                                                                  geneToPartners);
        // Do another method
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(clusterGenes,
                                                                        nodeToEdges);
        Map<String, Integer> pairToDistance2 = new HashMap<String, Integer>();
        for (String pair : pairToPath.keySet()) {
            List<String> path = pairToPath.get(pair);
            pairToDistance2.put(pair, path.size() - 1);
        }
        // Want to compare two maps
        System.out.println("Size: " + pairToDistance1.size() + ", " + pairToDistance2.size());
        for (String pair : pairToDistance1.keySet()) {
            Integer dist1 = pairToDistance1.get(pair);
            // Want to get the second
            Integer dist2 = pairToDistance2.get(pair);
            if (dist2 == null) {
                int index = pair.indexOf("\t");
                String id1 = pair.substring(0, index);
                String id2 = pair.substring(index + 1);
                pair = id2 + "\t" + id1;
                dist2 = pairToDistance2.get(pair);
            }
            if (dist1 != dist2) {
                System.out.println(pair + " has wrong distance: " + dist1 + ", " + dist2);
            }
        }
    }
    
    /**
     * This method is used to calculate a minimum spanning tree based on hierarchical
     * clustering based on shortest path.
     * @throws IOException
     */
    @Test
    public void calculateMinimumSpanBasedOnCluster() throws IOException {
//      String clusterFileName = DIR_NAME + "ClusterListForNatureGBM.txt";
//      String cutoff = "2.87";
//      String outFiFileName = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical.txt";
//      String nodeAttFileName = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical.na";
//      String outSpanFileName = DIR_NAME + "SpanCluster17ForScienceGBMViaHierarchical.txt";
//      String clusterFileName = DIR_NAME + "ClusterList.txt";
//      String cutoff = "3.06";
//      String outFiFileName = DIR_NAME + "FIsInCluster28ForPancreaticViaHierarchical.txt";
//      String nodeAttFileName = DIR_NAME + "FIsInCluster28ForPancreaticViaHierarchical.na";
//      String outSpanFileName = DIR_NAME + "SpanCluster28ForPancreaticViaHierarchical.txt";
//        String clusterFileName = DIR_NAME + "ClusterListForSciencePancreaticUp111208.txt";
//        String cutoff = "3.15";
//        String outFiFileName = DIR_NAME + "FIsInCluster30ForSciencePancreaticViaHierarchicalUp111208.txt";
//        String nodeAttFileName = DIR_NAME + "FIsInCluster30ForSciencePancreaticiaHierarchicalUp111208.na";
//        String outSpanFileName = DIR_NAME + "SpanCluster30ForSciencePancreaticViaHierarchicalUp111208.txt";
        double[] cutoffs = new double[] {
//                3.077777777777778,
//                3.0625,
//                3.0961538461538463
//                1.0
//                2.95
//                2.5245098039215685,
//                3.2170542635658914
                3.4134615384615383,
                3.1936936936936937,
                3.1745283018867925
        };
        String clusterFileNames[] = new String[] {
                //                DIR_NAME + "ClusterListForNatureGBM042109.txt",
                //                DIR_NAME + "ClusterListForScienceGBM042109.txt",
                //                DIR_NAME + "ClusterListForSciencePancreatic042109.txt"
                //                DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt"
                //                OVARIAN_RESULT_DIR + "ClusterList051909.txt"
//                TCGA_GBM_DIR_NAME + "ClusterListForTCGAGBM100909.txt",
//                TCGA_GBM_DIR_NAME + "ClusterListForParsonsGBM100909.txt"
                BREAST_AND_COLORECTAL_DIR + "ClusterListForBreastCancer030210.txt",
                BREAST_AND_COLORECTAL_DIR + "ClusterListForColorectalCancer030210.txt",
                PANCREATIC_DIR + "ClusterList030210.txt"
        };
        String outFiFileNames[] = new String[] {
                //DIR_NAME + "FIsInClusterNautreGBMViaHierarchical042109.txt",
                //DIR_NAME + "FIsInClusterScienceGBMViaHierarchical042109.txt",
                //DIR_NAME + "FIsInClusterSciencePancreaticViaHierarchical042109.txt"
               // DIR_NAME + "FIsInClusterListForNatureGBMViaSPath051409.txt"
//              OVARIAN_RESULT_DIR + "FIsInClusterList051909.txt"
//                TCGA_GBM_DIR_NAME + "FIsInTCGAAlteredClusterViaHiCluster101209.txt",
//                TCGA_GBM_DIR_NAME + "FIsInParsonsAlteredClusterViaHiCluster101209.txt"
                BREAST_AND_COLORECTAL_DIR + "FIsInClusterListForBreastCancer030210_test.txt",
                BREAST_AND_COLORECTAL_DIR + "FIsInClusterListForColorectalCancer030210_test.txt",
                PANCREATIC_DIR + "FIsInClusterList030210_test.txt"
        };
        String outSpanFileNames[] = new String[] {
                //DIR_NAME + "SpanInClusterNautreGBMViaHierarchical042109.txt",
                //DIR_NAME + "SpanInClusterScienceGBMViaHierarchical042109.txt",
                //DIR_NAME + "SpanInClusterSciencePancreaticViaHierarchical042109.txt"
                //    DIR_NAME + "SpanInClusterListForNatureGBMViaSPath051409.txt"
                //    OVARIAN_RESULT_DIR + "SpanInClusterList051909.txt"
//                TCGA_GBM_DIR_NAME + "SpanInTCGAAlteredClusterViaHiCluster101209.txt",
//                TCGA_GBM_DIR_NAME + "SpanInParsonsAlteredClusterViaHiCluster101209.txt"
                BREAST_AND_COLORECTAL_DIR + "SpanInClusterListForBreastCancer030210_test.txt",
                BREAST_AND_COLORECTAL_DIR + "SpanInClusterListForColorectalCancer030210_test.txt",
                PANCREATIC_DIR + "SpanInClusterList030210_Test.txt"
        };
        String[] geneExpFileNames = new String[] {
                //TCGA_GBM_DIR_NAME + "GeneExp_t_value_051409.txt"
                //OVARIAN_DIR_NAME + "GeneExp_t_value.txt"
        };
        for (int i = 0; i < cutoffs.length; i++) {
            String clusterFileName = clusterFileNames[i];
            System.out.println(clusterFileName);
            String cutoff = cutoffs[i] + "";
            String outFiFileName = outFiFileNames[i];
            //String geneExpFileName = geneExpFileNames[i];
            String geneExpFileName = null;
            String outSpanFileName = outSpanFileNames[i];
            List<String> cluster = getClusterForAnalysis(clusterFileName, cutoff);
            calculateMinimumSpan(outFiFileName, 
                                 cluster,
                                 geneExpFileName);
        }
    }
    
    /**
     * This method is used to calculate minimum span for all genes without clustering first.
     * @throws IOException
     */
    @Test
    public void calculateMinimumSpanAmongGenes() throws IOException {
        List<List<String>> allLists = new ArrayList<List<String>>();
//        List<String> genes = loadNatureGBMGenesWithTwoOrMoreMutations();
//        allLists.add(genes);
//        genes = loadScienceGBMGenes();
//        allLists.add(genes);
//        genes = loadGenesWithTwoOrMoreMutations();
//        allLists.add(genes);
        List<String> genes = new ArrayList<String>(loadOvarianMutatedGenes());
        allLists.add(genes);
        String outFiFileNames[] = new String[] {
//                DIR_NAME + "FIsInNautreGBM042109.txt",
//                DIR_NAME + "FIsInScienceGBM042109.txt",
//                DIR_NAME + "FIsInSciencePancreatic042109.txt"
//                DIR_NAME + "FIsInNautreGBM051409.txt",
///                DIR_NAME + "FIsInScienceGBM050409.txt",
//                DIR_NAME + "FIsInSciencePancreatic050409.txt"
                OVARIAN_RESULT_DIR + "FIsInAllMutatedGenes051909.txt"
        };
        String outSpanFileNames[] = new String[] {
//                DIR_NAME + "SpanInNautreGBM051409.txt",
//                DIR_NAME + "SpanInScienceGBM050409.txt",
//                DIR_NAME + "SpanInSciencePancreatic051409.txt"
                OVARIAN_RESULT_DIR + "SpanInAllMutatedGenes051909.txt"
        };
        String[] geneExpFileNames = new String[] {
                //TCGA_GBM_DIR_NAME + "GeneExp_t_value_051409.txt"
                OVARIAN_DIR_NAME + "GeneExp_t_value.txt"
        };
        // Used to filter non graph genes
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        for (int i = 0; i < allLists.size(); i++) {
            String outFiFileName = outFiFileNames[i];
            //String outSpanFileName = outSpanFileNames[i];
            String geneExpFileName = geneExpFileNames[i];
            List<String> cluster = allLists.get(i);
            cluster.retainAll(fiGenes);
            System.out.println("Total genes: " + cluster.size());
            calculateMinimumSpan(outFiFileName, 
                                 cluster,
                                 geneExpFileName);
        }
    }
    
    protected void calculateMinimumSpan(String outFiFileName,
                                      String outSpanFileName,
                                      List<String> cluster) throws IOException {
        calculateMinimumSpan(outFiFileName,
                             cluster,
                             null);
    }
    

    private void calculateMinimumSpan(String outFiFileName,
                                      List<String> cluster,
                                      String geneExpFileName) throws IOException {
//        Map<String, Integer> pairToDistance = calculateDistances(cluster);
//        List<HierarchicalClusterNode> newClusters = hierarchicalCluster(cluster, 
//                                                            pairToDistance, 
//                                                            true);
//        int totalWeight = 0;
//        // Calculate total cluster weight
//        HierarchicalClusterNode newFirstNode = newClusters.get(0);
//        outputCluster(newFirstNode);
//        Set<HierarchicalClusterNode> allNodes = new HashSet<HierarchicalClusterNode>();
//        HierarchicalCluster hclust = new HierarchicalCluster();
//        hclust.grepAllClusters(newFirstNode, allNodes);
//        for (HierarchicalClusterNode node : allNodes)
//            totalWeight += node.pathDistance;
//        System.out.println("Total weight: " + totalWeight);
//        // These are used to calculate the actual path
//        Set<String> spanningFIs = generateSpanFromClusters(cluster,
//                                                           allNodes,
//                                                           geneExpFileName);
//        System.out.println("Total FIs: " + spanningFIs.size());
//        fu.setOutput(outSpanFileName);
//        for (String fi : spanningFIs)
//            fu.printLine(fi);
//        fu.close();
//        // Want to generate all fis
//        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
//        generateFIsForCluster(outFiFileName, 
//                              fis, 
//                              spanningFIs, 
//                              cluster);
        NetworkBuilderForGeneSet builder = new NetworkBuilderForGeneSet();
        builder.setFiFileName(GENE_FI_BIG_COMP_FILE_NAME);
        builder.setPathwayFiGeneFileName(R3Constants.RESULT_DIR + "FIsInGene_Pathway_041709.txt");
        Set<String> fisInNetwork = builder.constructFINetworkForGeneSet(cluster, geneExpFileName);
        fu.saveInteractions(fisInNetwork, outFiFileName);
    }

    public Map<String, Double> loadGeneExpTValue(String fileName) throws IOException {
        fu.setInput(fileName);
        Map<String, Double> rtn = new HashMap<String, Double>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            rtn.put(tokens[0],
                    new Double(tokens[1])); 
        }
        return rtn;
    }

    private Set<String> generateSpanFromClusters(Set<HierarchicalClusterNode> allClusters,
                                                 Map<String, List<String>> pairToPath) {
        Set<String> spanningFIs = new HashSet<String>();
        for (HierarchicalClusterNode node : allClusters) {
            if (node.getChildNode1() == null ||
                node.getChildNode2() == null)
                continue; // The bottom node
            List<String> shortestPath = calculateShortestPath(node.getChildNode1(),
                                                              node.getChildNode2(), 
                                                              pairToPath);
            List<String> fis1 = convertPathToFIs(shortestPath);
            spanningFIs.addAll(fis1);
        }
        return spanningFIs;
    }
    
    private List<String> convertPathToFIs(List<String> path) {
        List<String> fis = new ArrayList<String>();
        for (int i = 0; i < path.size() - 1; i ++) {
            String id1 = path.get(i);
            String id2 = path.get(i + 1);
            // Want to make sure all FIs in in prefixed order
            int compare = id1.compareTo(id2);
            if (compare < 0)
                fis.add(id1 + "\t" + id2);
            else
                fis.add(id2 + "\t" + id1);
        }
        return fis;
    }
    
    /**
     * This method is used to generate a node attribute file for the GBM genes.The shared genes
     * between two data sets are colored different from non-shared genes.
     * @throws IOException
     */
    @Test
    public void generateNodeAttributeFilesForGBM() throws IOException {
        List<String> sharedGenes = getSharedGBMGenesFromNatureAndScience();
        // parameters for a cancer
        String clusterFileName = DIR_NAME + "ClusterListForNatureGBM102908.txt";
        String cutoff = "3.18";
        // For node attributes
        String nodeAttFileName = DIR_NAME + "FIsInCluster57ForNatureGBMViaHierarchical102908.na";
//      String clusterFileName = DIR_NAME + "ClusterListForScienceGBM102908.txt";
//      String cutoff = "2.94";
//      // For node attributes
//      String nodeAttFileName = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical102908.na";        
        // Start processing
        Map<String, List<String>> headerToCluster = loadClusters(clusterFileName);
        List<String> cancerGenes = null;
        fu.setOutput(nodeAttFileName);
        fu.printLine("node.fillColor (class=java.lang.String)");
        for (String header : headerToCluster.keySet()) {
            if (!(header.contains(cutoff)))
                continue; // Want to sure the smallest cluster only
            cancerGenes = headerToCluster.get(header);
            for (String gene : cancerGenes) {
                if (sharedGenes.contains(gene))
                    fu.printLine(gene + " = 204, 255, 153");
                else
                    fu.printLine(gene + " = 0, 204, 204");
            }
        }
        fu.close();
    }
    
    @Test
    public void generateNodeAttributeFilesForHubs() throws IOException {
        // These statements are used for degree based permutation
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        Set<String> hubProteins = new HashSet<String>();
        Set<String> nonHubProteins = new HashSet<String>();
        sortProteinsForHubs(proteinToDegree, hubProteins, nonHubProteins);
        
        List<String> sharedGenes = getSharedGBMGenesFromNatureAndScience();
        // parameters for a cancer
        String clusterFileName = DIR_NAME + "ClusterListForNatureGBM102908.txt";
        String cutoff = "3.18";
        // For node attributes
        String nodeAttFileName = DIR_NAME + "FIsInCluster57ForNatureGBMForHubs.na";
//      String clusterFileName = DIR_NAME + "ClusterListForScienceGBM102908.txt";
//      String cutoff = "2.94";
//      // For node attributes
//      String nodeAttFileName = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical102908.na";        
        // Start processing
        List<String> cancerGenes = getClusterForAnalysis(clusterFileName, cutoff);
        fu.setOutput(nodeAttFileName);
        fu.printLine("node.fillColor (class=java.lang.String)");
        for (String gene : cancerGenes) {
            if (hubProteins.contains(gene))
                fu.printLine(gene + " = 153, 153, 255");
            else
                fu.printLine(gene + " = 0, 204, 204");
        }
        fu.close();
    }
    
    /**
     * This is a permutation test to check a random pick genes for hierarchical analysis.
     * @throws Exception
     */
    @Test
    public void randomTest() throws Exception {
        long time1 = System.currentTimeMillis();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> sampled = randomSampling(genes, 62);
        String outFileName = DIR_NAME + "PathForRandom62.txt";
        generateShortestPathTable(new ArrayList<String>(sampled), outFileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for randomg test: " + (time2 - time1));
    }
    
    /**
     * This method is used to check how many cancer genes are hubs, and how many non-hubs.
     * @throws Exception
     */
    @Test
    public void checkCancerGeneDistributionsBasedOnDegrees() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        // These statements are used for degree based permutation
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        List<Set<String>> bins = sortAllGenesBasedOnDegrees(proteinToDegree);
        System.out.println("Sizes of bins:");
        int totalSize = proteinToDegree.keySet().size();
        for (int i = 0; i < bins.size(); i++) {
            Set<String> bin = bins.get(i);
            System.out.printf("%d\t%2.3f%n",
                              bin.size(),
                              bin.size() / (double) totalSize);
        }
        double[] cutoffs = new double[] {
                3.18,
                2.94,
                3.15
        };
        String clusterFileNames[] = new String[] {
                DIR_NAME + "ClusterListForNatureGBM102908.txt",
                DIR_NAME + "ClusterListForScienceGBM102908.txt",
                DIR_NAME + "ClusterListForSciencePancreatic102908.txt"
        };
        for (int i = 0; i < cutoffs.length; i++) {
            // Need to get genes
            String clusterFileName = clusterFileNames[i];
            double cutoff = cutoffs[i];
            List<String> genes = getClusterForAnalysis(clusterFileName, cutoff + "");
            bins = sortAllGenesBasedOnDegreesInCancerGenes(proteinToDegree, genes);
            totalSize = genes.size();
            // Find how many in hubs, how many in non-hubs
            System.out.println(clusterFileName);
            List<Integer> sizeInBin = calculateNumbersInBins(genes, bins);
            for (int j = 0; j < sizeInBin.size(); j++) {
                System.out.printf("%d\t%2.3f%n",
                                  sizeInBin.get(j),
                                  sizeInBin.get(j) / (double) totalSize);
            }
        }
    }
    
    /**
     * This method is used to check the cancer gene distributions based on bins
     * generated based on cancer gene distributions. Actually this is a test method
     * to see if the generated bins are right.
     * @throws Exception
     */
    @Test
    public void checkCancerGeneDegrees() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        // These statements are used for degree based permutation
        Map<String, Integer> geneToDegree = InteractionUtilities.generateProteinToDegree(fis);
        double averageForAll = calculateAverageDegree(geneToDegree.keySet(), 
                                                      geneToDegree);
        System.out.println("Total average degree: " + averageForAll);
        int allGenesTotal = geneToDegree.size();
        double[] cutoffs = new double[] {
                2.5245098039215685,
                3.2170542635658914
//                1.0,
//                3.0625,
                //3.15
        };
        String clusterFileNames[] = new String[] {
                //                DIR_NAME + "ClusterListForNatureGBM102908.txt",
                //                DIR_NAME + "ClusterListForScienceGBM102908.txt",
                //                DIR_NAME + "ClusterListForSciencePancreatic102908.txt"
                //                DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt",
                //                DIR_NAME + "ClusterListForScienceGBM042109.txt"
                TCGA_GBM_DIR_NAME + "ClusterListForTCGAGBM100909.txt",
                TCGA_GBM_DIR_NAME + "ClusterListForParsonsGBM100909.txt"
        };
        for (int i = 0; i < cutoffs.length; i++) {
            // Need to get genes
            String clusterFileName = clusterFileNames[i];
            double cutoff = cutoffs[i];
            List<String> genes = getClusterForAnalysis(clusterFileName, cutoff + "");
            double average = calculateAverageDegree(genes, geneToDegree);
            double pvalue = ttestForDegree(averageForAll, genes, geneToDegree);
            List<Set<String>> bins = sortAllGenesBasedOnDegreesInCancerGenes(geneToDegree,
                                                                             genes);
            int totalSize = genes.size();
            // Find how many in hubs, how many in non-hubs
            System.out.println(clusterFileName);
            System.out.println("Average degree: " + average + "; pvalue: " + pvalue);
            List<Integer> sizeInBin = calculateNumbersInBins(genes, bins);
            for (int j = 0; j < sizeInBin.size(); j++) {
                Set<String> bin = bins.get(j);
                System.out.printf("%d\t%2.3f\t%d\t%2.3f%n",
                                  bin.size(),
                                  bin.size() / (double) allGenesTotal,
                                  sizeInBin.get(j),
                                  sizeInBin.get(j) / (double) totalSize);
            }
        }
    }
    
    private double ttestForDegree(double mu,
                                  Collection<String> genes,
                                  Map<String, Integer> geneToDegree) throws MathException {
        double[] values = new double[genes.size()];
        int index = 0;
        for (String gene : genes) {
            Integer degree = geneToDegree.get(gene);
            values[index] = degree;
            index ++;
        }
        return TestUtils.tTest(mu, values);
    }
    
    private double calculateAverageDegree(Collection<String> genes,
                                          Map<String, Integer> geneToDegree) {
        int total = 0;
        for (String gene : genes) {
            Integer degree = geneToDegree.get(gene);
            total += degree;
        }
        return (double) total / genes.size();
    }
    
    /**
     * This test will target for hierarchical clustering but with average shortest path calcualted
     * out. The results have not been sorted.
     * @throws Exception
     */
    @Test
    public void permutationTestForHierarchicalClusterShortestBasedOnDegree() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // These statements are used for degree based permutation
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        // Want to permutated pre-selected GBM genes
        //Set<String> allGBMGenes = loadAllNatureGBM601Genes();
        //proteinToDegree.keySet().retainAll(allGBMGenes);
        System.out.println("Total sampling genes: " + proteinToDegree.size());
        //List<Set<String>> bins = sortAllGenesBasedOnDegrees(proteinToDegree);
        //Set<String> hubProteins = new HashSet<String>();
        //Set<String> nonHubProteins = new HashSet<String>();
        //sortProteinsForHubs(proteinToDegree, hubProteins, nonHubProteins);
        // For hierarchical clustering
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        int permutationNumber = 1000;
        // Load Science GBM
        //List<String> genes = loadScienceGBMGenes();
        // Parsons altered genes
        NatureGBMAnalyzer gbmAnalyzer = new NatureGBMAnalyzer();
        Map<String, Set<String>> sampleToAlteredGenes = gbmAnalyzer.getScienceGBMSampleToAlteredGenes();
        List<String> genes = new ArrayList<String>(gbmAnalyzer.selectGenesInSamples(2, sampleToAlteredGenes));
        genes.retainAll(allGenes);
        System.out.println("Science GBM: " + genes.size());
        String outFileName = DIR_NAME + "PermutationTestForClusterSciGBMBasedOnDegreeSelfBinWithAvgDist_1000_010410.txt";
        double cutoff = 3.2170542635658914; 
        List<Set<String>> bins = sortAllGenesBasedOnDegreesInCancerGenes(proteinToDegree,
                                                                         genes);
        List<Integer> sizeInBins = calculateNumbersInBins(genes, bins);
        System.out.println("Starting permutation...");
        StringBuilder builder = new StringBuilder();
        for (int j = 0; j < permutationNumber; j ++) {
            System.out.println("Permutation " + j + "...");
            //long time1 = System.currentTimeMillis();
            //          Set<String> sampled = new HashSet<String>();
            //          Set<String> sampled1 = randomSampling(hubProteins, hubInGenes.size());
            //          Set<String> sampled2 = randomSampling(nonHubProteins, nonHubInGenes.size());
            //          sampled.addAll(sampled1);
            //          sampled.addAll(sampled2);
            Set<String> sampled = samplingBasedOnBins(sizeInBins, bins);
            List<String> sampledList = new ArrayList<String>(sampled);
            Map<String, Integer> pairToDistance = calculateDistances(sampledList, 
                                                                     bfs, 
                                                                     geneToPartners);
            List<HierarchicalClusterNode> clusters = hierarchicalCluster(sampled, 
                                                                         pairToDistance,
                                                                         false);
            HierarchicalClusterNode root = clusters.get(0);
            HierarchicalClusterNode largestGoodNode = searchLargestNode(root, cutoff);
            double ratio = (double) largestGoodNode.ids.size() / root.ids.size();
            double path = calculateShortestPath(sampledList,
                                                bfs, 
                                                nodeToEdges);
            builder.append((j + 1) + "\t" + ratio + "\t" + path);
            builder.append("\n");
        }
        fu.setOutput(outFileName);
        fu.printLine(builder.toString());
        fu.close();
    }
    
    @Test
    public void permutationTestForAvgShortestPathBasedOnDegree() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        // These statements are used for degree based permutation
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        // Want to permutated pre-selected GBM genes
//        Set<String> allGBMGenes = loadAllNatureGBM601Genes();
//        proteinToDegree.keySet().retainAll(allGBMGenes);
//        System.out.println("Total sampling genes: " + proteinToDegree.size());
         
//        List<Set<String>> bins = sortAllGenesBasedOnDegrees(proteinToDegree);
//        Set<String> hubProteins = new HashSet<String>();
//        Set<String> nonHubProteins = new HashSet<String>();
//        sortProteinsForHubs(proteinToDegree, hubProteins, nonHubProteins);
        // For hierarchical clustering
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        int permutationNumber = 100;
        double[] cutoffs = new double[] {
                //3.18,
                //2.94,
                //3.15
                //1.0
                //2.95
                //3.0625
                //2.5606060606060606
//                2.5245098039215685,
                3.2170542635658914
//                3.4134615384615383,
//                3.473684210526316,
//                3.1936936936936937,
//                3.03125,
//                3.1745283018867925
        };
        String clusterFileNames[] = new String[] {
                //DIR_NAME + "ClusterListForNatureGBM102908.txt",
                //DIR_NAME + "ClusterListForScienceGBM102908.txt",
//              TCGA_GBM_DIR_NAME + "ClusterListForTCGAGBM100909.txt",
              TCGA_GBM_DIR_NAME + "ClusterListForParsonsGBM100909.txt"
                //DIR_NAME + "ClusterListForSciencePancreatic102908.txt"
                //DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt"
                //OVARIAN_RESULT_DIR + "ClusterList051909.txt"
//                BREAST_AND_COLORECTAL_DIR + "ClusterListForBreastCancer030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "ClusterListForBreastCancer030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "ClusterListForColorectalCancer030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "ClusterListForColorectalCancer030210.txt",
//                PANCREATIC_DIR + "ClusterList030210.txt"
        };
        // Permutation test based on bins sortAllGenesBasedOnDegrees(Map<String, Integer>)
//        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeBin_NatureGBM_1000_111408.txt",
//                DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeBin_ScienceGBM_1000_111408.txt",
//                DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeBin_SciencePancreatic_1000_111408.txt"
//        };
        // Permutation test based on bins from sortAllGenesBasedOnDegreesInCancerGenes
        // (Map<String, Integer>, List<String>)
        String[] outFileNames = new String[] {
//                BREAST_AND_COLORECTAL_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_Breast_All_Samples_1000_030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_Breast_Disc_Samples_1000_030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_Colorectal_All_Samples_1000_030210.txt",
//                BREAST_AND_COLORECTAL_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_Colorectal_Disc_Samples_1000_030210.txt",
//                PANCREATIC_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_Pancreatic_1000_030210.txt"
//                TCGA_GBM_DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_TCGAGBM_100_042810.txt",
                TCGA_GBM_DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_ParsonsGBM_100_042810.txt",
                //DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_ScienceGBM_100_093008.txt",
//                DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_ScienceGBM_1000_100809.txt",
//                DIR_NAME + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_SciencePancreatic_1000_111508.txt"
//                DIR_NAME + "PermutationTestForAvgShortestPathBaseadOnDegreeSelfBinPre_NatureGBMVisaSPath_Cut_1_081009.txt"
//                OVARIAN_RESULT_DIR + "PermutationTestForAvgShortestPathBasedOnDegreeSelfBin_1000_051909.txt",
        };
        for (int i = 0; i < cutoffs.length; i++) {
            // Need to get genes
            String clusterFileName = clusterFileNames[i];
            double cutoff = cutoffs[i];
            List<String> genes = getClusterForAnalysis(clusterFileName, 
                                                       cutoff + "");
            // Just in case if there is anything is not correct.
            genes.retainAll(proteinToDegree.keySet());
            List<Set<String>> bins = sortAllGenesBasedOnDegreesInCancerGenes(proteinToDegree, 
                                                                             genes);
            List<Integer> numbersInBins = calculateNumbersInBins(genes, 
                                                                 bins);
//            // Find how many in hubs, how many in non-hubs
//            List<String> hubInGenes = new ArrayList<String>(genes);
//            hubInGenes.retainAll(hubProteins);
//            List<String> nonHubInGenes = new ArrayList<String>(genes);
//            nonHubInGenes.removeAll(hubInGenes);
            System.out.println("Starting permutation...");
            List<Double> values = new ArrayList<Double>();
            for (int j = 0; j < permutationNumber; j ++) {
                long time1 = System.currentTimeMillis();
//                Set<String> sampled = new HashSet<String>();
//                Set<String> sampled1 = randomSampling(hubProteins, hubInGenes.size());
//                Set<String> sampled2 = randomSampling(nonHubProteins, nonHubInGenes.size());
//                sampled.addAll(sampled1);
//                sampled.addAll(sampled2);
                Set<String> sampled = samplingBasedOnBins(numbersInBins, 
                                                          bins);
                List<String> sampledList = new ArrayList<String>(sampled);
                double value = calculateShortestPath(new ArrayList<String>(sampled),
                                                     bfs, 
                                                     nodeToEdges);
                values.add(value);
                long time2 = System.currentTimeMillis();
                System.out.println(j + ": " + (time2 - time1));
            }
            Collections.sort(values);
            String outFileName = outFileNames[i];
            fu.setOutput(outFileName);
            for (double ratio : values)
                fu.printLine(ratio + "");
            statDescRandomValues(values, fu);
            fu.close();
        }
    }
    
    private Set<String> samplingBasedOnBins(List<Integer> sizeInBins,
                                            List<Set<String>> bins) {
        Set<String> sample = new HashSet<String>();
        for (int i = 0; i < sizeInBins.size(); i++) {
            int size = sizeInBins.get(i);
            Set<String> bin = bins.get(i);
            Set<String> sampled = randomSampling(bin, size);
            sample.addAll(sampled);
        }
        return sample;
    }
    
    /**
     * Permutation test for average shortest path based on a pre-defined list.
     * @throws Exception
     */
    @Test
    public void permutationTestForAverageShortestPathOnPreList() throws Exception {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        // Get the test genes
        Set<String> preListGenes = loadAllNatureGBM601Genes();
        System.out.println("Total Pre-list genes: " + preListGenes.size());
        preListGenes.retainAll(proteins);
        System.out.println("Total pre-list genes in the network: " + preListGenes.size());
        int permutation = 1000;
        int[] sampleSizes = new int[] {
                50,
//                30, 
//                17
//                21
        };
        String[] fileNames = new String[] {
                DIR_NAME + "PermutationTestForAverageShortestPathPreList_50_1000_081009.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_30_1000_103008.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_17_1000_103008.txt"
//                DIR_NAME + "PermutationTestForAverageShortestPath_21_1000_110308.txt"
        };
        for (int i = 0; i < sampleSizes.length; i++) {
            int sampleSize = sampleSizes[i];
            String fileName = fileNames[i];
            List<Double> values = new ArrayList<Double>();
            for (int j = 0; j < permutation; j++) {
                Set<String> samples = randomSampling(preListGenes, sampleSize);
                double value = calculateShortestPath(new ArrayList<String>(samples),
                                                     bfs, 
                                                     nodeToEdges);
                values.add(value);
                //System.out.println(i + ": " + value);
            }
            Collections.sort(values);
            fu.setOutput(fileName);
            for (Double value : values)
                fu.printLine(value + "");
            fu.close();
            //for (Double value : values)
            //    System.out.println(value);
        }
        System.out.println("Done: permutationTestForAverageShortestPathOnPreList()");
    }
    
    /**
     * This method is used to do permutation test for average shortest for a fixed number of gene set.
     * @throws Exception
     */
    @Test
    public void permutationTestForAverageShortestPath() throws Exception {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> proteinList = new ArrayList<String>(proteins);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        int permutation = 100;
        int[] sampleSizes = new int[] {
//                57,
//                30, 
//                17
//                21
//                11
//                18
//                82,
//                20,
//                80,
//                61
                46, 
                71
        };
        String[] fileNames = new String[] {
//                DIR_NAME + "PermutationTestForAverageShortestPath_57_1000_103008.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_30_1000_103008.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_17_1000_103008.txt"
//                DIR_NAME + "PermutationTestForAverageShortestPath_21_1000_110308.txt"
//                R3Constants.RESULT_DIR + "PermutationTestForAverageShortestPath_11_1000_041709.txt"
//                R3Constants.RESULT_DIR + "PermutationTestForAverageShortestPath_18_1000_080709.txt"
                DIR_NAME + "PermutationTestForAverageShortestPath_46_100_042810.txt",
                DIR_NAME + "PermutationTestForAverageShortestPath_71_100_042810.txt"
//                DIR_NAME + "PermutationTestForAverageShortestPath_82_1000_030210.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_20_1000_030210.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_80_1000_030210.txt",
//                DIR_NAME + "PermutationTestForAverageShortestPath_61_1000_030210.txt"
        };
        for (int i = 0; i < sampleSizes.length; i++) {
            int sampleSize = sampleSizes[i];
            String fileName = fileNames[i];
            List<Double> values = new ArrayList<Double>();
            for (int j = 0; j < permutation; j++) {
                Set<String> samples = randomSampling(proteinList, sampleSize);
                double value = calculateShortestPath(new ArrayList<String>(samples),
                                                     bfs, 
                                                     nodeToEdges);
                values.add(value);
                System.out.println(j + ": " + value);
            }
            Collections.sort(values);
            fu.setOutput(fileName);
            for (Double value : values)
                fu.printLine(value + "");
            statDescRandomValues(values, fu);
            fu.close();
            //for (Double value : values)
            //    System.out.println(value);
        }
    }

    private void statDescRandomValues(List<Double> values,
                                      FileUtility fu) throws IOException {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (Double v : values)
            stat.addValue(v);
        fu.printLine("Min: " + stat.getMin());
        fu.printLine("Max: " + stat.getMax());
        fu.printLine("Mean: " + stat.getMean());
        fu.printLine("SD: " + stat.getStandardDeviation());
    }
    
    private void sortProteinsForHubs(Map<String, Integer> proteinToDegree,
                                     Set<String> hubs,
                                     Set<String> nonHubs) {
        for (String protein : proteinToDegree.keySet()) {
            Integer degree = proteinToDegree.get(protein);
            if (degree > HUB_CUTOFF)
                hubs.add(protein);
            else
                nonHubs.add(protein);
        }
    }
    
    /**
     * This helper method is used to create bins based on degree distributions in the cancer genes.
     * @param proteinToDegree
     * @param cancerGenes
     * @return
     */
    private List<Set<String>> sortAllGenesBasedOnDegreesInCancerGenes(Map<String, Integer> proteinToDegree,
                                                                     List<String> cancerGenes) {
        List<Set<String>> bins = new ArrayList<Set<String>>();
        // Get a list of degrees for the cancer genes
        Set<Integer> degreesInCancerGenes = new HashSet<Integer>();
        for (String gene : cancerGenes) {
            Integer degree = proteinToDegree.get(gene);
            degreesInCancerGenes.add(degree);
        }
        List<Integer> degreeList = new ArrayList<Integer>(degreesInCancerGenes);
        // Do a reverse sorting
        Collections.sort(degreeList, new Comparator<Integer>() {
            public int compare(Integer degree1, Integer degree2) {
                return degree2 - degree1;
            };
        });
        for (int i = 0; i < degreeList.size(); i++) {
            bins.add(new HashSet<String>());
        }
//         Check the degrees
        System.out.println("Check degrees:");
        for (Integer degree : degreeList)
            System.out.println(degree);
        // Sort genes based on degrees
        for (String gene : proteinToDegree.keySet()) {
            Integer degree = proteinToDegree.get(gene);
            // Find which bin this gene should be
            for (int i = 0; i < degreeList.size(); i++) {
                Integer cutoff = degreeList.get(i);
                if (degree >= cutoff) {
                    Set<String> bin = bins.get(i);
                    bin.add(gene);
                    break;
                }
            }
        }
        return bins;
    }
    
    /**
     * A helper method to create protein bins based on protein degrees.
     * @param proteinToDegree
     * @return
     */
    private List<Set<String>> sortAllGenesBasedOnDegrees(Map<String, Integer> proteinToDegree) {
        List<Set<String>> bins = new ArrayList<Set<String>>();
        // Manually created bins
        List<int[]> thresholds = new ArrayList<int[]>();
        thresholds.add(new int[]{666, 220}); // 0-1% (approximately)
        thresholds.add(new int[]{219, 176}); // 1-2%
        thresholds.add(new int[]{175, 152}); // 2-3%
        thresholds.add(new int[]{151, 136}); // 3-4%
        thresholds.add(new int[]{135, 128}); // 4-5%
        thresholds.add(new int[]{127, 105}); // 5-7%
        thresholds.add(new int[]{104, 87}); // 7-9%
        thresholds.add(new int[]{86, 74}); // 9-11%
        thresholds.add(new int[]{73, 36}); // 11-20%
        thresholds.add(new int[]{35, 19}); // 20-30%
        thresholds.add(new int[]{18, 11}); // 30-40%
        thresholds.add(new int[]{10, 7}); // 40-50%
        thresholds.add(new int[]{6, 3}); // 50-70%
        thresholds.add(new int[]{2, 1}); // 70-100%
        for (int i = 0; i < thresholds.size(); i++) {
            bins.add(new HashSet<String>());
        }
        for (String protein : proteinToDegree.keySet()) {
            int degree = proteinToDegree.get(protein);
            // Find which bin this protein should be
            for (int i = 0; i < thresholds.size(); i++) {
                int[] values = thresholds.get(i);
                if (degree <= values[0] && degree >= values[1]) {
                    Set<String> bin = bins.get(i);
                    bin.add(protein);
                    break;
                }
            }
        }
        // Just a quick check
//        for (Set<String> bin : bins)
//            System.out.println(bin.size());
        return bins;
    }
    
    private List<Integer> calculateNumbersInBins(List<String> genes,
                                                 List<Set<String>> bins) {
        List<Integer> sizes = new ArrayList<Integer>(bins.size());
        for (int i = 0; i < bins.size(); i++) {
            Set<String> bin = bins.get(i);
            Set<String> copy = new HashSet<String>(genes);
            copy.retainAll(bin);
            sizes.add(copy.size());
        }
        return sizes;
    }
    
    @Test
    public void permutationTestForHierarchicalClusterBasedOnDegree() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // These statements are used for degree based permutation
        Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(fis);
        // Want to permutated pre-selected GBM genes
//        Set<String> allGBMGenes = loadAllNatureGBM601Genes();
//        proteinToDegree.keySet().retainAll(allGBMGenes);
        System.out.println("Total sampling genes: " + proteinToDegree.size());
        // For hierarchical clustering
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        int permutationNumber = 1000;
        // Load Science Pancreatic 
        List<List<String>> allLists = new ArrayList<List<String>>();
        // Ovarian genes
//        List<String> genes = new ArrayList<String>(loadOvarianMutatedGenes());
//        genes.retainAll(allGenes);
//        allLists.add(genes);
////        List<String> genes = loadGenesWithTwoOrMoreMutations();
////        genes.retainAll(allGenes);
////        System.out.println("Science Pancreatic genes: " + genes.size());
////        allLists.add(genes);
        // Load Nature GBM
//        List<String> genes = loadNatureGBMGenesWithTwoOrMoreMutations();
//        genes.retainAll(allGenes);
//        System.out.println("Nature GBM: " + genes.size());
//        allLists.add(genes);
//        // Load Science GBM
//        List<String> genes = loadScienceGBMGenes();
//        genes.retainAll(allGenes);
//        System.out.println("Science GBM: " + genes.size());
//        allLists.add(genes);
//        // Permutation test based on bins generated by method, sortAllGenesBasedOnDegrees(Map<String, Integer>).
//        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForClusterSciPanBasedOnDegreeBin_1000_111408.txt",
//                DIR_NAME + "PermutationTestForClusterNatGBMBasedOnDegreeBin_1000_111408.txt",
//                DIR_NAME + "PermutationTestForClusterSciGBMBasedOnDegreeBin_1000_111408.txt"
//        };
        // TCGA altered genes
        NatureGBMAnalyzer gbmAnalyzer = new NatureGBMAnalyzer();
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Map<String, Set<String>> sampleToAlteredGenes = gbmAnalyzer.getSampleToAlteredGenes(samples);
        List<String> genes = new ArrayList<String>(gbmAnalyzer.selectGenesInSamples(2, sampleToAlteredGenes));
        genes.retainAll(proteinToDegree.keySet());
        allLists.add(genes);
        // Parsons altered genes
        sampleToAlteredGenes = gbmAnalyzer.getScienceGBMSampleToAlteredGenes();
        genes = new ArrayList<String>(gbmAnalyzer.selectGenesInSamples(2, sampleToAlteredGenes));
        genes.retainAll(proteinToDegree.keySet());
        allLists.add(genes);
        // Permutation test based on bins generated by method, 
        // sortAllGenesBasedOnDegreesInCancerGenes()
        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForClusterSciPanBasedOnDegreeSelfBin_1000_111508.txt",
//                DIR_NAME + "PermutationTestForClusterNatGBMBasedOnDegreeSelfBin_1000_111508.txt",
//                DIR_NAME + "PermutationTestForClusterSciGBMBasedOnDegreeSelfBin_1000_111508.txt"
//              DIR_NAME + "PermutationTestForClusterViaSPathSciPanBasedOnDegreeSelfBin_1000_042109.txt",
//              DIR_NAME + "PermutationTestForClusterViaSPathNatGBMBasedOnDegreeSelfBinPre_1000_cutoff_1_081009.txt",
//              DIR_NAME + "PermutationTestForClusterViaSPathSciGBMBasedOnDegreeSelfBin_1000_042109.txt"
//                OVARIAN_RESULT_DIR + "PermutationTestForClusterBasedOnDegreeSelfBin_1000_051909.txt",
//                DIR_NAME + "PermutationTestForClusterSciGBMBasedOnDegreeSelfBin_1000_081009.txt"
                TCGA_GBM_DIR_NAME + "PermutationTestForClusteringSelfBin_252_98_1000_101209.txt",
                TCGA_GBM_DIR_NAME + "PermutationTestForClusteringSelfBin_322_65_1000_101209.txt",
        };
        double[] cutoffs = new double[] {
                //2.0,
                //1.0,
                //2.0
                //2.95
                //3.0625
                2.5245098039215685,
                3.2170542635658914
        };
        for (int i = 0; i < allLists.size(); i++) {
            // Need to get genes
            genes = allLists.get(i);
            List<Set<String>> bins = sortAllGenesBasedOnDegreesInCancerGenes(proteinToDegree,
                                                                             genes);
            List<Integer> sizeInBins = calculateNumbersInBins(genes, bins);
            System.out.println("Starting permutation...");
            List<Double> values = new ArrayList<Double>();
            double cutoff = cutoffs[i];
            for (int j = 0; j < permutationNumber; j ++) {
                Set<String> sampled = samplingBasedOnBins(sizeInBins, bins);
                List<String> sampledList = new ArrayList<String>(sampled);
                Map<String, Integer> pairToDistance = calculateDistances(sampledList, 
                                                                         bfs, 
                                                                         geneToPartners);
                List<HierarchicalClusterNode> clusters = hierarchicalCluster(sampled, 
                                                                             pairToDistance,
                                                                             false);
                HierarchicalClusterNode root = clusters.get(0);
                HierarchicalClusterNode largestGoodNode = searchLargestNode(root, cutoff);
                double ratio = (double) largestGoodNode.ids.size() / root.ids.size();
                values.add(ratio);
                //System.out.println(j + ": " + ratio);
                //long time2 = System.currentTimeMillis();
                //System.out.println(i + ": " + (time2 - time1));
            }
            // Do a reverse sorting
            Collections.sort(values, new Comparator<Double>() {
                public int compare(Double value1, Double value2) {
                    return value2.compareTo(value1);
                }
            });
            String outFileName = outFileNames[i];
            fu.setOutput(outFileName);
            for (double ratio : values)
                fu.printLine(ratio + "");
            statDescRandomValues(values, fu);
            fu.close();
        }
        System.out.println("Done: permutationTestForHierarchicalClusterBasedOnDegree()");
    }
    
    /**
     * Permutation test for hierarchilca cluster based on a pre-list.
     * @throws Exception
     */
    @Test
    public void permutationTestForHierarchicalClusterBasedOnPreList() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        // Load prelist
        Set<String> preList = loadAllNatureGBM601Genes();
        preList.retainAll(genes);
        
        int permutationNumber = 1000;
        double[] cutoffs = new double[] {
                1.0,
               // 2.94,
               // 3.15
        };
        int[] sampleSizes = new int[] {
                68,
               // 25,
               // 44
        };
        String[] outFileNames = new String[] {
                //DIR_NAME + "PermutationTestForClusteringPreList_318_68_1000_112108_1.txt",
                DIR_NAME + "PermutationTestForClusteringPreListViaSPath_1_68_1000_081009.txt",
               // DIR_NAME + "PermutationTestForClustering_294_25_1000_103008.txt",
               // DIR_NAME + "PermutationTestForClustering_315_44_1000_103008.txt"
        };
        for (int j = 0; j < cutoffs.length; j++) {
            double cutoff = cutoffs[j];
            int sampleSize = sampleSizes[j];
            String outFileName = outFileNames[j];
            System.out.println("Starting permutation...");
            List<Double> values = new ArrayList<Double>();
            for (int i = 0; i < permutationNumber; i ++) {
                //long time1 = System.currentTimeMillis();
                Set<String> sampled = randomSampling(preList, sampleSize);
                List<String> sampledList = new ArrayList<String>(sampled);
                Map<String, Integer> pairToDistance = calculateDistances(sampledList, 
                                                                         bfs, 
                                                                         geneToPartners);
                List<HierarchicalClusterNode> clusters = hierarchicalCluster(sampled, 
                                                                 pairToDistance,
                                                                 false);
                HierarchicalClusterNode root = clusters.get(0);
                HierarchicalClusterNode largestGoodNode = searchLargestNode(root, cutoff);
                double ratio = (double) largestGoodNode.ids.size() / root.ids.size();
                values.add(ratio);
                //long time2 = System.currentTimeMillis();
                //System.out.println(i + ": " + (time2 - time1));
            }
            Collections.sort(values);
            fu.setOutput(outFileName);
            for (double ratio : values)
                fu.printLine(ratio + "");
            fu.close();
        }
        System.out.println("Done: permutationTestForHierarchicalClusterBasedOnPreList()");
    }
    
    /**
     * This method is used to do a multiple permutation test.
     * @throws Exception
     */
    @Test
    public void permutationTestForHierarchicalCluster() throws Exception {
        // Set up fixed variables
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        int permutationNumber = 100;
        double[] cutoffs = new double[] {
                //                3.18,
                //                2.94,
                //                3.15
                //                3.0625
                //                1.0
                //                3.25
//                2.5245098039215685,
                3.2170542635658914
//                3.4134615384615383,
//                3.473684210526316
//                3.1936936936936937,
//                3.03125
//                3.1745283018867925
        };
        int[] sampleSizes = new int[] {
//                98,
                65,
//                25,
//                44
//                113,
//                28
//                113,
//                25
//                79
        };
        String[] outFileNames = new String[] {
//                DIR_NAME + "PermutationTestForClustering_317_79_1000_030210.txt",
//              DIR_NAME + "PermutationTestForClustering_319_113_1000_030210.txt",
//              DIR_NAME + "PermutationTestForClustering_303_25_1000_030210.txt",
//                DIR_NAME + "PermutationTestForClustering_341_113_1000_030210.txt",
//                DIR_NAME + "PermutationTestForClustering_348_28_1000_030210.txt",
//                DIR_NAME + "PermutationTestForClustering_252_98_100_042810.txt",
                DIR_NAME + "PermutationTestForClustering_322_65_100_042810.txt",
//                DIR_NAME + "PermutationTestForClustering_294_25_1000_103008.txt",
//                DIR_NAME + "PermutationTestForClustering_315_44_1000_103008.txt"
//                DIR_NAME + "PermutationTestForClustering_30625_25_1000_081009.txt"
//                DIR_NAME + "PermutationTestForClusteringViaSPath_1_68_1000_081009.txt"
        };
        for (int j = 0; j < cutoffs.length; j++) {
            double cutoff = cutoffs[j];
            int sampleSize = sampleSizes[j];
            String outFileName = outFileNames[j];
            System.out.println("Starting permutation...");
            List<Double> values = new ArrayList<Double>();
            for (int i = 0; i < permutationNumber; i ++) {
                //long time1 = System.currentTimeMillis();
                Set<String> sampled = randomSampling(genes, sampleSize);
                List<String> sampledList = new ArrayList<String>(sampled);
                Map<String, Integer> pairToDistance = calculateDistances(sampledList, 
                                                                         bfs, 
                                                                         geneToPartners);
                List<HierarchicalClusterNode> clusters = hierarchicalCluster(sampled, 
                                                                 pairToDistance,
                                                                 false);
                HierarchicalClusterNode root = clusters.get(0);
                HierarchicalClusterNode largestGoodNode = searchLargestNode(root, cutoff);
                double ratio = (double) largestGoodNode.ids.size() / root.ids.size();
                values.add(ratio);
                //long time2 = System.currentTimeMillis();
                //System.out.println(i + ": " + (time2 - time1));
            }
            Collections.sort(values);
            fu.setOutput(outFileName);
            for (double ratio : values)
                fu.printLine(ratio + "");
            statDescRandomValues(values, fu);
            fu.close();
        }
    }
    
    protected HierarchicalClusterNode searchLargestNode(HierarchicalClusterNode root,
                                          double cutoff) {
        List<HierarchicalClusterNode> goodNodes = new ArrayList<HierarchicalClusterNode>();
        searchGoodNodes(root,
                        cutoff,
                        goodNodes);
        // Need to find the largest nodes
        HierarchicalClusterNode largestGoodNode = goodNodes.get(0);
        for (HierarchicalClusterNode node : goodNodes) {
            if (largestGoodNode.ids.size() < node.ids.size())
                largestGoodNode = node;
        }
        return largestGoodNode;
    }
    
    private void searchGoodNodes(HierarchicalClusterNode checkNode,
                                 double cutoff,
                                 List<HierarchicalClusterNode> goodNode) {
        if (checkNode.pathDistance <= cutoff) {
            goodNode.add(checkNode);
            return;
        }
        // Need to check left and right nodes
        searchGoodNodes(checkNode.getChildNode1(), cutoff, goodNode);
        searchGoodNodes(checkNode.getChildNode2(), cutoff, goodNode);
    }
    
    protected Set<String> randomSampling(Collection<String> genes,
                                       int size) {
        return MathUtilities.randomSampling(genes, size);
    }
    
    /**
     * This method is used to check the coverage of the list of mutated genes.
     * @throws IOException
     */
    @Test
    public void checkGeneCoverage() throws Exception {
        // This is based on gene names
        Set<String> fis = fu.loadInteractions(GENE_FI_FILE_NAME);
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        //Collection<String> mutatedGenes = loadGenesWithTwoOrMoreMutations();
        //Collection<String> mutatedGenes = loadNatureGBMGenesWithTwoOrMoreMutations();
        //Collection<String> mutatedGenes = loadScienceGBMGenes();
        Collection<String> mutatedGenes = loadOvarianMutatedGenes();
        int total = mutatedGenes.size();
        mutatedGenes.retainAll(genes);
        int hits = mutatedGenes.size();
        double ratio = (double) hits/total;
        System.out.println("Based on gene names:");
        System.out.println("Total: " + total);
        System.out.println("In the FI network: " + hits);
        System.out.println("ratio: " + ratio);
//        // This is based on UniProt ids
//        fis = fu.loadInteractions(R3Constants.INTERACTION_FILE_NAME);
//        genes = InteractionUtilities.grepIDsFromInteractions(fis);
//        mutatedGenes = loadIdsWithTwoOrMoreMutations();
//        total = mutatedGenes.size();
//        mutatedGenes.retainAll(genes);
//        hits = mutatedGenes.size();
//        ratio = (double) hits/total;
//        System.out.println("Based on UniProt Ids:");
//        System.out.println("Total: " + total);
//        System.out.println("In the FI network: " + hits);
//        System.out.println("ratio: " + ratio);
    }
    
    /**
     * This method is used to do a hierarchical clustering.
     * Note on Aug 7, 2009: A completely new implementation to remove an intermediate file.
     * @throws Exception
     */
    @Test
    public void hierarchicalClusterMutatedGenes() throws Exception {
        // Genes to be clustering
        List<String> alteredGenes = loadScienceGBMGenes();
        hierarchicalClusterGenes(alteredGenes);
    }

    protected void hierarchicalClusterGenes(List<String> alteredGenes) throws IOException {
        // Loading Fis
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total FI genes in the biggest network: " + fiGenes.size());
        System.out.println("Total altered genes: " + alteredGenes.size());
        alteredGenes.retainAll(fiGenes);
        System.out.println("Genes in FI network: " + alteredGenes.size());
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        Map<String, Integer> pairToDistance = calculateDistances(alteredGenes, 
                                                                 bfs, 
                                                                 geneToPartners);
        List<HierarchicalClusterNode> clusters = hierarchicalCluster(alteredGenes, 
                                                                     pairToDistance, 
                                                                     false);
        HierarchicalClusterNode cluster = clusters.get(0);
        outputCluster(cluster);
    }
    
    protected void outputCluster(HierarchicalClusterNode root) {
        System.out.println("hierarchical layouting...");
        new HierarchicalCluster().outputCluster("", root);
    }
    
    private double distance(String id1,
                            String id2,
                            Map<String, Double> pairToPath) {
        String key = id1 + "\t" + id2;
        Double path = pairToPath.get(key);
        if (path == null) {
            key = id2 + "\t" + id1;
            path = pairToPath.get(key);
        }
        return path;
    }

    public List<HierarchicalClusterNode> hierarchicalCluster(Collection<String> ids,
                                                                Map<String, ? extends Number> pairToPath, // Use String - Integer or String - double
                                                                boolean useShortestPath) {
        // Make sure pairToPath should be String - Double
        final Map<String, Double> pairToPath1 = new HashMap<String, Double>();
        for (Object key : pairToPath.keySet()) {
            Object value = pairToPath.get(key);
            pairToPath1.put(key.toString(),
                            new Double(value.toString()));
        }
        DistanceCalculator distanceCalculator = new DistanceCalculator() {
            
            public double calculateDistance(String id1, String id2) {
                return distance(id1, id2, pairToPath1);
            }
        };
        
        HierarchicalCluster cluster = new HierarchicalCluster();
        cluster.setDistanceCalculator(distanceCalculator);
        if (useShortestPath)
            cluster.setMethod(HierarchicalCluster.ClusterDistanceMethod.SINGLE);
        else
            cluster.setMethod(HierarchicalCluster.ClusterDistanceMethod.AVERAGE);
        HierarchicalClusterNode top = cluster.cluster(ids);
        List<HierarchicalClusterNode> rtn = new ArrayList<HierarchicalClusterNode>();
        rtn.add(top);
        return rtn;
    }
    
    @Test
    public void checkConnectionInClusters() throws Exception {
        Graph<String, DefaultEdge> graph = new GraphAnalyzer().createGraph(R3Constants.RESULT_DIR + "FI73InGene_061008.txt");
        String clusterFileName = DIR_NAME + "ClusterList.txt";
        Map<String, List<String>> headerToGenes = loadClusters(clusterFileName);
        Set<String> mutatedGenes = getSomaticMutatedGenes();
        for (String header : headerToGenes.keySet()) {
            List<String> genes = headerToGenes.get(header);
            System.out.println(header);
            for (int i = 0; i < genes.size() - 1; i++) {
                String gene1 = genes.get(i);
                for (int j = i + 1; j < genes.size(); j++) {
                    String gene2 = genes.get(j);
                    List<DefaultEdge> path = DijkstraShortestPath.findPathBetween(graph, gene1, gene2);
                    Set<String> genesInPaht = getGenesFromPath(path);
                    if (!mutatedGenes.containsAll(genesInPaht))
                        System.out.println(gene1 + " " + gene2 + " is not connected by mutated genes");
                }
            }
        }
    }
    
    private Set<String> getGenesFromPath(List<DefaultEdge> path) {
        Set<String> genes = new HashSet<String>();
        for (DefaultEdge edge : path) {
            String tmp = edge.toString();
            String sub = tmp.substring(1, tmp.length() - 1);
            String[] tokens = sub.split(":");
            genes.add(tokens[0].trim());
            genes.add(tokens[1].trim());
        }
        return genes;
    }
    
    private Set<String> getSomaticMutatedGenes() throws IOException {
        //String fileName = DIR_NAME + "SomaticMutationsInDiscoveryScreen.txt";
        String fileName = DIR_NAME + "SomaticMutationsInPrevalenceScreen.txt";
        fu.setInput(fileName);
        Set<String> genes = new HashSet<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            String type = tokens[6];
            if (!type.equals("Synonymous"))
                genes.add(gene);
        }
        fu.close();
        return genes;
    }
    
    private Map<String, Set<String>> loadSampleToGenes() throws IOException {
        String[] fileNames = new String[] {
                DIR_NAME + "SomaticMutationsInDiscoveryScreen.txt",
                DIR_NAME + "SomaticMutationsInPrevalenceScreen.txt"
        };
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            String line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                String gene = tokens[0];
                // Used to use lower case since the problem in the discover screen file
                String sample = tokens[2].toLowerCase();
                String type = tokens[6];
                if (type.equals("Synonymous"))
                    continue;
                Set<String> set = sampleToGenes.get(sample);
                if (set == null) {
                    set = new HashSet<String>();
                    sampleToGenes.put(sample, set);
                }
                set.add(gene);
            }
            fu.close();
        }
        return sampleToGenes;
    }
    
    /**
     * This method is used to check the distributions of genes in the cancer samples.
     * @throws IOException
     */
    @Test
    public void checkGeneInSampleDistribution() throws IOException {
        Map<String, Set<String>> sampleToGenes = loadSampleToGenes();
        final Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        List<String> genes = new ArrayList<String>(geneToSamples.keySet());
        Collections.sort(genes, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Set<String> samples1 = geneToSamples.get(gene1);
                Set<String> samples2 = geneToSamples.get(gene2);
                return samples2.size() - samples1.size();
            }
        });
        for (String gene : genes) {
            Set<String> samples = geneToSamples.get(gene);
            System.out.println(gene + "\t" + samples.size());
        }
        System.out.println("Total samples: " + sampleToGenes.size());
    }
    
    /**
     * Load all 601 genes for TCGA cancer re-sequencing project.
     * @return
     * @throws IOException
     */
    protected Set<String> loadAllNatureGBM601Genes() throws IOException {
        String fileName = DIR_NAME + "TCGA601Genes.txt";
        Set<String> genes = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[1]);
        }
        fu.close();
        return genes;
    }
    
    private Set<String> loadNatureGlioblastomaGenes() throws IOException {
        String fileName = DIR_NAME + "GlioblastomaMutationTable.txt";
        Set<String> genes = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[0]);
        }
        fu.close();
        return genes;
    }
    
    @Test
    public void checkSamplesInMutationTable() throws IOException {
        String fileName = TCGA_GBM_DIR_NAME + "TCGA_GBM_Level3_Somatic_Mutations_08.28.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> samples = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            samples.add(tokens[15]);
        }
        System.out.println("Total samples: " + samples.size());
    }
    
    protected List<String> loadNatureGBMGenesWithTwoOrMoreMutations() throws IOException {
        Map<String, Integer> geneToCount = loadNatureGBMMutationTable();
        List<String> rtn = new ArrayList<String>();
        for (String gene : geneToCount.keySet()) {
            Integer count = geneToCount.get(gene);
            if (count > 1)
                rtn.add(gene);
        }
        return rtn;
    }

    /**
     * Load Nature GBM mutation table.
     * @return key: gene names; value: number of mutation.
     * @throws IOException
     */
    public Map<String, Integer> loadNatureGBMMutationTable()
            throws IOException {
        String fileName = DIR_NAME + "GlioblastomaMutationTable.txt";
        Map<String, Integer> geneToCount = new HashMap<String, Integer>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Integer count = geneToCount.get(tokens[0]);
            if (count == null)
                geneToCount.put(tokens[0], 1);
            else
                geneToCount.put(tokens[0], ++count);
        }
        fu.close();
        return geneToCount;
    }
    
    public Set<String> loadOvarianMutatedGenes() throws IOException {
        //String fileName = OVARIAN_DIR_NAME + "Ovarian_phaseI_somatic_mutation_table.050809.txt";
        String fileName = OVARIAN_DIR_NAME + "Ovarian_phaseI_somatic_mutation_table.051309.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> mutatedGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            mutatedGenes.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total mutated genes: " + mutatedGenes.size());
        return mutatedGenes;
    }
    
    /**
     * This method is used to check ids in the shortest path file.;
     * @throws Exception
     */
    @Test
    public void checkIdsInPathFile() throws Exception {
        String fileName = DIR_NAME + "PathAmongIdsWithTwoOrMore.txt";
        fu.setInput(fileName);
        String line = null;
        Set<String> ids = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            ids.add(tokens[0]);
            ids.add(tokens[1]);
        }
        fu.close();
        System.out.println("Total ids: " + ids.size());
    }
    
    /**
     * This method is used to generate a pair-wise shortest path among 
     * a list of genes.
     * @throws Exception
     */
    @Test
    public void generateShortestDistanceTable() throws Exception {
        List<List<String>> allLists = new ArrayList<List<String>>();
        //List<String> genes = loadGenesWithTwoOrMoreMutations();
        //allLists.add(genes);
        //List<String> genes = loadLungCancerGenes();
        //List<String> genes = loadNatureGBMGenesWithTwoOrMoreMutations();
        List<String> genes = new ArrayList<String>(loadOvarianMutatedGenes());
        //allLists.add(genes);
        //genes = loadScienceGBMGenes();
        //List<String> genes1 = loadScienceGBMGenes();
        //genes.addAll(genes1);
        allLists.add(genes);
        String[] fileNames = new String[] {
//                DIR_NAME + "PathAmongSciencePancreaticUp111208.txt",
//                DIR_NAME + "PathwayAmongNatureGBMUp111208.txt",
//                DIR_NAME + "PathAmongScienceGBMGenesUp111208.txt"
//              DIR_NAME + "PathAmongSciencePancreatic042109.txt",
//              DIR_NAME + "PathwayAmongNatureGBM042109.txt",
//              DIR_NAME + "PathAmongScienceGBMGenes042109.txt"   
//                DIR_NAME + "PathAmongBothGBMGenes042109.txt"
                OVARIAN_RESULT_DIR + "PathAmongMutatedGenes051909.txt"
        };
        for (int i = 0; i < allLists.size(); i++) {
            genes = allLists.get(i);
            System.out.println("Total genes: " + genes.size());
            String outFileName = fileNames[i];
            generateShortestPathTable(genes, outFileName);
        }
    }
    
    /**
     * Use this method to calculate shortest path from files generated from method
     * generateShortestDistanceTable().
     * @throws IOException
     */
    @Test
    public void calculateShortestPathAmongGenes() throws IOException {
        String[] fileNames = new String[] {
//              DIR_NAME + "PathAmongSciencePancreaticUp111208.txt",
//              DIR_NAME + "PathwayAmongNatureGBMUp111208.txt",
//              DIR_NAME + "PathAmongScienceGBMGenesUp111208.txt"
//                DIR_NAME + "PathAmongSciencePancreatic042109.txt",
//                DIR_NAME + "PathwayAmongNatureGBM042109.txt",
//                DIR_NAME + "PathAmongScienceGBMGenes042109.txt"   
//                DIR_NAME + "PathAmongBothGBMGenes042109.txt"
                OVARIAN_RESULT_DIR + "PathAmongMutatedGenes051909.txt"
        };
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            System.out.println(fileName);
            String line = null;
            int total = 0;
            int count = 0;
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                total += Integer.parseInt(tokens[2]);
                count ++;
            }
            System.out.println("Average path: " + (double)total/count);
        }
    }
    
    /**
     * This method is used to geneate FIs among all shortest paths for ids
     * in a specified cluster.
     * @throws Exception
     */
    @Test
    public void generateFIsBasedOnShortestPath() throws Exception {
        String clusterFileName = DIR_NAME + "ClusterList.txt";
        String cutoff = "3.06";
        String outFileName = DIR_NAME + "FIsFromAllShortestPancreatic.txt";
        List<String> genes = getClusterForAnalysis(clusterFileName, 
                                                   cutoff);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(genes, nodeToEdges);
        Set<String> fisInCluster = new HashSet<String>();
        for (String pair : pairToPath.keySet()) {
            List<String> path = pairToPath.get(pair);
            for (int i = 0; i < path.size() - 1; i++) {
                String name1 = path.get(i);
                String name2 = path.get(i + 1);
                int compare = name1.compareTo(name2);
                if (compare < 0)
                    fisInCluster.add(name1 + "\t" + name2);
                else if (compare > 0)
                    fisInCluster.add(name2 + "\t" + name1);
            }
        }
        fu.setOutput(outFileName);
        for (String fi : fisInCluster)
            fu.printLine(fi);
        fu.close();
    }

    private void generateShortestPathTable(List<String> genes,
                                           String outFileName) throws IOException {
        Map<String, Integer> pairToDistance = calculateDistances(genes);
        // Call the following statement after the first line since fu is used
        // in the first method call.
        fu.setOutput(outFileName);
        for (String pair : pairToDistance.keySet()) {
            int length = pairToDistance.get(pair);
            System.out.println(pair + "\t" + length);
            if (length == Integer.MAX_VALUE)
                fu.printLine(pair + "\tNA");
            else
                fu.printLine(pair + "\t" + length);
        }
        fu.close();
    }
    
    private Map<String, Integer> calculateDistances(List<String> genes) throws IOException {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        return calculateDistances(genes, bfs, geneToPartners);
    }

    /**
     * Calcualte average distance for a list of genes.
     * @param genes
     * @param bfs
     * @param geneToPartners
     * @return
     */
    public Map<String, Integer> calculateDistances(List<String> genes,
                                                      BreadthFirstSearch bfs,
                                                      Map<String, Set<String>> geneToPartners) {
        Map<String, Integer> pairToDistance = new HashMap<String, Integer>();
        List<String> targets = new ArrayList<String>();
        for (int i = 0; i < genes.size() - 1; i++) {
            String id1 = genes.get(i);
            if (!(geneToPartners.containsKey(id1)))
                continue;
            for (int j = i + 1; j < genes.size(); j++) {
                String id2 = genes.get(j);
                if (!(geneToPartners.containsKey(id2)))
                    continue;
                targets.add(id2);
                //int length = bfs.getDistance(id1, id2, geneToPartners);
                //pairToDistance.put(id1 + "\t" + id2, length);
            }
            Map<String, Integer> targetToDist = bfs.getDistances(id1,
                                                                 targets,
                                                                 geneToPartners);
            for (String target : targetToDist.keySet()) {
                Integer dist = targetToDist.get(target);
                pairToDistance.put(id1 + "\t" + target, dist);
            }
            targets.clear();
        }
        return pairToDistance;
    }
    
    /**
     * This method is used to annotate a set of genes.
     * @throws Exception
     */
    @Test
    public void annotateGeneSets() throws Exception {
        //List<String> names = loadScienceGBMGenes();
        //List<String> names = loadNatureGBMGenesWithTwoOrMoreMutations();
        //Collection<String> names = loadOvarianMutatedGenes();
        String geneList = "PSG1, VNN1, TNFRSF10C, NOVA2, ZNF800, " +
        		"NIPBL, BBC3, HNRNPUL2, NUCB2, LIPN, ERGIC2, MLLT4, " +
        		"FBXL20, CHAD, LOC641367, SLC4A8, FEN1, ZNF426, TARS, C12orf44";
        String[] tmp = geneList.split(", ");
        List<String> names = Arrays.asList(tmp);
        annotateNames(names);
        //Map<String, Set<String>> nameToIds = getGeneNameToIdsMap();
        //annotateIds(names, nameToIds);
    }

    private void annotateIds(Collection<String> names,
                             Map<String, Set<String>> nameToIds) throws Exception {
        // Need to convert gene names to gene ids
        Set<String> ids = new HashSet<String>();
        for (String name : names) {
            Set<String> tmp = nameToIds.get(name);
            if (tmp != null)
                ids.addAll(tmp);
        }
        PathwayBasedAnnotator topicAnalyzer = new PathwayBasedAnnotator();
//        topicAnalyzer.setHopNumber(0);
        topicAnalyzer.setPValueThreshold(1.0);
        topicAnalyzer.annotateGenesWithFDR(ids, AnnotationType.Pathway);
        //System.out.println();
        //System.out.println("Annotation using clustered pathways:");
        //topicAnalyzer.clusteredTopicAnnotate(ids, System.out);
    }
    
    protected void annotateNames(Collection<String> names) throws Exception {
        PathwayBasedAnnotator topicAnalyzer = new PathwayBasedAnnotator();
        //topicAnalyzer.setPValueThreshold(0.01);
        topicAnalyzer.setFDRThreshold(0.25);
        System.out.println("Annotate with the original pathways:");
        topicAnalyzer.annotateGenesWithFDR(names, AnnotationType.Pathway);
//        System.out.println("\nAnnotate with GO BP:");
//        topicAnalyzer.annotateGenesUsingGOWithFDR(names, "BP");
//        System.out.println("\nAnnotate with GO CC:");
//        topicAnalyzer.annotateGenesUsingGOWithFDR(names, "CC");
//        System.out.println("\nAnnotate with GO MF:");
//        topicAnalyzer.annotateGenesUsingGOWithFDR(names, "MF");
//        System.out.println();
//        System.out.println("Annotate with the clustered pathways:");
//        topicAnalyzer.topicAnnotateWithNames(names, true, System.out);
    }
                               
    /**
     * This method is used to annotate a mini-span sub graph generated from method 
     * calculateMinimumSpanBasedOnCluster().
     * @throws Exception
     */
    @Test
    public void annotateSpanSubGraph() throws Exception {
        String[] fileNames = new String[] {
                //DIR_NAME + "SpanCluster30ForSciencePancreaticViaHierarchical102908.txt",
                //DIR_NAME + "SpanCluster17ForScienceGBMViaHierarchical102908.txt",
                //DIR_NAME + "SpanCluster57ForNatureGBMViaHierarchical102908.txt"
                //DIR_NAME + "FIsInClusterScienceGBMViaHierarchical042109.txt"
                TCGA_GBM_DIR_NAME + "FIsInTCGAAlteredClusterViaHiCluster101209.txt",
                TCGA_GBM_DIR_NAME + "FIsInParsonsAlteredClusterViaHiCluster101209.txt"
        };
        for (String fileName : fileNames) {
            System.out.println("File: " + fileName);
            Set<String> fis = fu.loadInteractions(fileName);
            Set<String> names = InteractionUtilities.grepIDsFromInteractions(fis);
            annotateNames(names);
            System.out.println();
        }
    }
    
    /**
     * This method is used to annotate FIs extracted from pathways.
     * @throws Exception
     */
    @Test
    public void annotateFIs() throws Exception {
//        String input = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical102908.txt";
//        String output = DIR_NAME + "FIsInCluster17ForScienceGBMViaHierarchical102908WithAnnot.txt";
//        String input = DIR_NAME + "FIsInCluster57ForNatureGBMViaHierarchical102908.txt";
//        String output = DIR_NAME + "FIsInCluster57ForNatureGBMViaHierarchical102908WithAnnot.txt";
//        String input = DIR_NAME + "SpanCluster57ForNatureGBMViaHierarchical102908.txt";
//        String output = DIR_NAME + "SpanCluster57ForNatureGBMViaHierarchical102908WithAnnot.txt";
//        String input = DIR_NAME + "FIsInCluster30ForSciencePancreaticViaHierarchical102908.txt";
//        String output = DIR_NAME + "FIsInCluster30ForSciencePancreaticViaHierarchical102908Annot.txt";
        String input = DIR_NAME + "FIsInNautreGBM042109.txt";
        String output = DIR_NAME + "FIsInNautreGBM042109Annot.txt";
        FISourceTypeReader typeQuery = new FISourceTypeReader();
        Session session = typeQuery.initSession().openSession();
        Set<String> fis = fu.loadInteractions(input);
        fu.setOutput(output);
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index);
            String name2 = fi.substring(index + 1);
            Map<String, Set<String>> pairToTypes = typeQuery.queryTypes(name1, 
                                                                        name2, 
                                                                        session);
//            if (pairToTypes.size() > 1) {
//                for (String pair : pairToTypes.keySet()) {
//                    System.out.println(pair + ": " + generateType(pairToTypes.get(pair)));
//                }
//                throw new IllegalStateException("Two directions for types: " + fi);
//            }
            // If there are two directions, pick one based on annotated information
//            if (pairToTypes.size() > 1) {
//                for (Iterator<String> it = pairToTypes.keySet().iterator(); it.hasNext();) {
//                    String pair = it.next();
//                    Set<String> types = pairToTypes.get(pair);
//                    String type = generateType(types);
//                    System.out.println(type);
//                    if (!type.contains("activate") && 
//                        !type.contains("catalyze") &&
//                        !(type.contains("inhibit") && !type.contains("inhibition"))) {
//                        pairToTypes.remove(pair);
//                        System.out.println(" Removed!");
//                        break; // An arbitrary way to keep one direction
//                    }
//                }
//            }
            for (String pair : pairToTypes.keySet()) {
                Set<String> types = pairToTypes.get(pair);
                String type = typeQuery.generateType(types);
                if (pair.equals(fi))
                    fu.printLine(name1 + "\t" + type + "\t" + name2);
                else
                    fu.printLine(name2 + "\t" + type + "\t" + name1);
            }
        }
        fu.close();
        session.close();
    }
    
    /**
     * Used to annotate a list of clusters in a file.
     * @throws Exception
     */
    @Test
    public void annotateClusters() throws Exception {
//        String clusterFileName = DIR_NAME + "ClusterListForSciencePancreatic102908.txt";
//        String cutoff = "3.15";
        String clusterFileName = TCGA_GBM_DIR_NAME + "ClusterListForParsonsGBM100909.txt";
        String cutoff = "3.2170542635658914";
//        String clusterFileName = DIR_NAME + "ClusterListForNatureGBMViaSPath042109.txt";
//        String cutoff = "1.0";
//        String clusterFileName = OVARIAN_RESULT_DIR + "ClusterList051909.txt";
//        String cutoff = "2.95";
        Map<String, List<String>> headerToIds = loadClusters(clusterFileName);
        for (String tmp : headerToIds.keySet()) {
            if (!tmp.contains(cutoff))
                continue;
            System.out.println(tmp);
            List<String> names = headerToIds.get(tmp);
            annotateNames(names);
        }
    }

    private Map<String, Set<String>> getGeneNameToIdsMap() throws Exception {
        HibernateFIReader hibernateAnalyzer = new HibernateFIReader();
        Map<String, Set<String>> nameToIds = hibernateAnalyzer.generateProteinNameToAccession();
        return nameToIds;
    }
    
    /**
     * This method is used to calculate average shortest path among ids in clusters.
     * @throws Exception
     */
    @Test
    public void calculateAverageShortestPathForClusters() throws Exception {
        String[] fileNames = new String[] {
                //"ClusterListForScienceGBM.txt",
                //"ClusterList.txt",
                //"ClusterListForGlioblastoma.txt",
//                "ClusterListForSciencePancreatic102908.txt",
//                "ClusterListForScienceGBM102908.txt",
//                "ClusterListForNatureGBM102908.txt"
//                "ClusterListForSciencePancreaticUp111208.txt",
//                "ClusterListForScienceGBMUp111208.txt",
//                "ClusterListForNatureGBMUp111208.txt"
//                "ClusterListForScienceGBMViaSPath042109.txt",
//                "ClusterListForNatureGBMViaSPath042109.txt",
//                "ClusterListForSciencePancreaticViaSPath042109.txt",
//                "ClusterListForSciencePancreatic042109.txt",
//                "ClusterListForNatureGBM042109.txt",
                "ClusterListForTCGAGBM100909.txt",
                //"ClusterListForScienceGBM042109.txt"
                "ClusterListForParsonsGBM100909.txt"
//                "ClusterList051909.txt"
        };
        for (String clusterFileName : fileNames) {
            System.out.println(clusterFileName);
            clusterFileName = TCGA_GBM_DIR_NAME + clusterFileName;
            calculateAverageShortestPathInClusters(clusterFileName);
        }
    }

    protected void calculateAverageShortestPathInClusters(String clusterFileName)
            throws IOException {
        Map<String, List<String>> headerToIds = loadClusters(clusterFileName);
        //Map<String, List<String>> headerToIds = loadClusters(OVARIAN_RESULT_DIR + clusterFileName);
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(GENE_FI_BIG_COMP_FILE_NAME);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        long time1 = System.currentTimeMillis();
        for (String header : headerToIds.keySet()) {
            List<String> ids = headerToIds.get(header);
            //Set<String> allGBMGenes = loadAllNatureGBM601Genes();
            //ids.retainAll(allGBMGenes);
            double value = calculateShortestPath(ids, bfs, nodeToEdges);
            System.out.println(header);
            System.out.println(value);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time used: " + (time2 - time1));
        System.out.println();
    }
    
    /**
     * Calculate average shortest path for a set of genes
     * @param ids
     * @param bfs
     * @param nodeToEdges
     * @return
     */
    protected double calculateShortestPath(List<String> ids,
                                           BreadthFirstSearch bfs,
                                           Map<TreeNode, List<Edge>> nodeToEdges) {
        Map<String, List<String>> pairToPath = bfs.generateShortestPath(ids,
                                                                        nodeToEdges);
        int total = 0;
        for (String pair : pairToPath.keySet()) {
            total += pairToPath.get(pair).size() - 1;
        }
        double value = (double) total / pairToPath.size();
        return value;
    }
    
    /**
     * This method is used to check shared genes in clusters from two different cancers.
     * @throws IOException
     */
    @Test
    public void checkSharedClusterGenesFromDiffCancers() throws IOException {
        Map<String, List<String>> pancreaticClusters = loadClusters(DIR_NAME + "ClusterList.txt");
        Map<String, List<String>> glioblastomaClusters = loadClusters(DIR_NAME + "ClusterListForGlioblastoma.txt");
        List<String> pancreaticCluster = null;
        for (String header : pancreaticClusters.keySet()) {
            if (header.contains("3.06")) {
                pancreaticCluster = pancreaticClusters.get(header);
            }
        }
        List<String> glioblastomaCluster = null;
        for (String header : glioblastomaClusters.keySet()) {
            if (header.contains("3.06"))
                glioblastomaCluster = glioblastomaClusters.get(header);
        }
        // Check shared genes
        glioblastomaCluster.retainAll(pancreaticCluster);
        System.out.println("Total shared: " + glioblastomaCluster.size());
        for (String gene : glioblastomaCluster)
            System.out.println(gene);
    }
    
    /**
     * This method is used to check genes in the clusters.
     * @throws IOException
     */
    @Test
    public void checkGenesInClusters() throws IOException {
        String clusterFileName = DIR_NAME + "ClusterList.txt";
        Map<String, List<String>> headerToCluster = loadClusters(clusterFileName);
        Set<String> mutatedGenes = getSomaticMutatedGenes();
        // Check how many mutated genes in the clusters
        for (String header : headerToCluster.keySet()) {
            System.out.println(header);
            List<String> genes = headerToCluster.get(header);
            for (String gene : genes)
                System.out.println(gene);
            Set<String> copy = new HashSet<String>(mutatedGenes);
            int total = copy.size();
            copy.retainAll(genes);
            int overlap = copy.size();
            double ratio = (double) overlap / total;
            System.out.println(overlap + " / " + total + " (" + ratio + ")");
        }
    }
    
    /**
     * This method is used to check the sample coverage in the found clusters
     * @throws IOException
     */
    @Test
    public void checkSamplesInClusters() throws IOException {
        String clusterFileName = DIR_NAME + "ClusterList.txt";
        Map<String, List<String>> headerToCluster = loadClusters(clusterFileName);
        Map<String, Set<String>> sampleToGenes = loadSampleToGenes();
//        List<String> samples = new ArrayList<String>(sampleToGenes.keySet());
//        Collections.sort(samples);
//        for (String sample : samples)
//            System.out.println(sample);
        for (String header : headerToCluster.keySet()) {
            List<String> genes = headerToCluster.get(header);
            List<String> coveredSamples = new ArrayList<String>();
            for (String sample : sampleToGenes.keySet()) {
                Set<String> sampleGenes = sampleToGenes.get(sample);
                sampleGenes.retainAll(genes);
                if (sampleGenes.size() > 0)
                    coveredSamples.add(sample);
            }
            System.out.println(header);
            double ratio = (double) coveredSamples.size() / sampleToGenes.size();
            System.out.println(coveredSamples.size() + " / " + sampleToGenes.size() + " (" + ratio + ")");
        }
    }

    protected Map<String, List<String>> loadClusters(String fileName) throws IOException {
        Map<String, List<String>> headerToIds = new HashMap<String, List<String>>();
        fu.setInput(fileName);
        String line = null;
        String header = null;
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                continue; // Empty line
            if (line.startsWith("#")) {
                header = line;
            }
            else {
                String[] tokens = line.split(", ");
                List<String> ids = new ArrayList<String>();
                for (String token : tokens)
                    ids.add(token);
                headerToIds.put(header, ids);
            }
        }
        fu.close();
        return headerToIds;
    }
    
    
    @Test
    public void checkGeneToUniProtMap() throws Exception {
        List<String> genes = loadNatureGBMGenesWithTwoOrMoreMutations();
        //HibernateFIAnalyzer fiAnalyzer = new HibernateFIAnalyzer();
        //Map<String, Set<String>> nameToIds = fiAnalyzer.generateProteinNameToAccession();
        UCSCDataAnalyzer analyzer = new UCSCDataAnalyzer();
        Map<String, String> geneToUni = analyzer.mapToUniProt(genes);
        for (String gene : geneToUni.keySet()) {
            String uni = geneToUni.get(gene);
            System.out.println(gene + "->" + uni);
        }
    }
    
    private List<String> loadScienceGenes(List<String> genes, String fileName)
            throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[0]);
        }
        fu.close();
        return genes;
    }
    
    /**
     * This method is used to load all genes from two mutation tables (table 3
     * and table 4). Genes contained in CNV are not included.
     * @return
     * @throws IOException
     */
    protected Set<String> loadAllScienceGBMMutatedGenes() throws IOException {
        String gbmDir = DIR_NAME + "GBMDataset/";
        String[] fileNames = new String[] {
                gbmDir + "1164382tableS3.txt",
                gbmDir + "1164382tableS4.txt"
        };
        Set<String> genes = new HashSet<String>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            // There are two annotation lines
            String line = fu.readLine();
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                genes.add(tokens[0]);
            }
            fu.close();
        }
        return genes;
    }
    
    protected Map<String, Set<String>> loadAllScienceSampleToMutatedGenes() throws IOException {
        String gbmDir = DIR_NAME + "GBMDataset/";
        String[] fileNames = new String[] {
                gbmDir + "1164382tableS3.txt",
                gbmDir + "1164382tableS4.txt"
        };
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            // There are two annotation lines
            String line = fu.readLine();
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                String gene = tokens[0];
                String sample = tokens[2];
                String mutationType = tokens[6];
                if (mutationType.equals("Synonymous"))
                    continue; // Don't want to include these genes
                Set<String> genes = sampleToGenes.get(sample);
                if (genes == null) {
                    genes = new HashSet<String>();
                    sampleToGenes.put(sample, genes);
                }
                genes.add(gene);
            }
            fu.close();
        }
        return sampleToGenes;
    }
    
    /**
     * This method will return altered genes in both CNV and mutation tables. Note: the samples
     * in the returned Map are for 22 samplese in the discovery screen since only these 22 samples
     * have CNV data set available.
     * @return
     * @throws IOException
     */
    protected Map<String, Set<String>> getScienceGBMSampleToAlteredGenes() throws IOException {
        Map<String, Set<String>> sampleToMutatedGenes = loadAllScienceSampleToMutatedGenes();
        Map<String, Set<String>> sampleToCNVGenes = loadAllScienceSampleToCNVGenes();
        Set<String> samples = new HashSet<String>();
        samples.addAll(sampleToCNVGenes.keySet());
        Map<String, Set<String>> merged = new HashMap<String, Set<String>>();
        for (String sample : samples) {
            Set<String> genes = new HashSet<String>();
            Set<String> mutated = sampleToMutatedGenes.get(sample);
            if (mutated != null)
                genes.addAll(mutated);
            Set<String> cnvGenes = sampleToCNVGenes.get(sample);
            if (cnvGenes != null)
                genes.addAll(cnvGenes);
            merged.put(sample, genes);
        }
        return merged;
    }
    
    protected Map<String, Set<String>> loadAllScienceSampleToCNVGenesViaUCSC() throws IOException {
        String gbmDir = DIR_NAME + "GBMDataset/";
        String[] fileNames = new String[] {
                "1164382tableS6.txt",
                "1164382tableS5.txt"
        };  
        Map<String, Set<String>> sampleToCNVGenes = new HashMap<String, Set<String>>();
        Set<String> samples = new HashSet<String>();
        Map<String, List<String[]>> geneToCoordinates = new UCSCDataAnalyzer().loadGeneNameToCoordinates();
        for (String fileName : fileNames) {
            fu.setInput(gbmDir + fileName);
            String line = fu.readLine();
            fu.readLine();
            while ((line = fu.readLine()) != null) {
          //      System.out.println(line);
                String[] tokens = line.split("\t");
                String sample = tokens[0];
                samples.add(sample);
                Set<String> genes = sampleToCNVGenes.get(sample);
                if (genes == null) {
                    genes = new HashSet<String>();
                    sampleToCNVGenes.put(sample, genes);
                }
                // Want to find genes
                String chr = "chr" + tokens[1];
                String token = tokens[2].substring(1, tokens[2].length() - 1);
                token = token.replaceAll(",", "");
                int start = Integer.parseInt(token);
                token = tokens[3].substring(1, tokens[3].length() - 1);
                token = token.replaceAll(",", "");
                int end = Integer.parseInt(token);
                for (String gene : geneToCoordinates.keySet()) {
                    List<String[]> list = geneToCoordinates.get(gene);
                    for (String[] coord : list) {
                        if (!coord[0].equals(chr))
                            continue;
                        int start0 = Integer.parseInt(coord[1]);
                        int end0 = Integer.parseInt(coord[2]);
                        // Check if there any overlap
                        if (end > start0 & start < end0) {
                            genes.add(gene);
                        }
                    }
                }
            }
            fu.close();
        }
        System.out.println("Total samples in CNVs: " + samples.size());
        for (String sample : sampleToCNVGenes.keySet()) {
            Set<String> genes = sampleToCNVGenes.get(sample);
            System.out.println(sample + ": " + genes.size());
        }
        return sampleToCNVGenes;

    }
    
    protected Map<String, Set<String>> loadAllScienceSampleToCNVGenes() throws IOException {
        String gbmDir = DIR_NAME + "GBMDataset/";
        String[] fileNames = new String[] {
                gbmDir + "1164382tableS6.txt",
                gbmDir + "1164382tableS5.txt"
        };  
        return loadAllScienceSampleToCNVGenes(fileNames, 4, 6);
    }

    protected Map<String, Set<String>> loadAllScienceSampleToCNVGenes(String[] fileNames,
                                                                      int geneCol,
                                                                      int otherGeneCol) throws IOException {
        Map<String, Set<String>> sampleToCNVGenes = new HashMap<String, Set<String>>();
        Set<String> samples = new HashSet<String>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            String line = fu.readLine();
            fu.readLine();
            while ((line = fu.readLine()) != null) {
                if (line.trim().length() == 0)
                    break;
//                System.out.println(line);
                String[] tokens = line.split("\t");
                String sample = tokens[0];
                samples.add(sample);
                String gene = tokens[geneCol];
                Set<String> genes = sampleToCNVGenes.get(sample);
                if (genes == null) {
                    genes = new HashSet<String>();
                    sampleToCNVGenes.put(sample, genes);
                }
                if (!gene.equals("Unknown")) {
                    genes.add(gene);
                }
                // Check other genes
                if (tokens.length == (otherGeneCol + 1) && tokens[otherGeneCol].length() > 0) {
                    // Do some simple parsing
                    String token = null;
                    if (tokens[otherGeneCol].startsWith("\""))
                        token = tokens[otherGeneCol].substring(1, tokens[otherGeneCol].length() - 1); // Remove quotation
                    else
                        token = tokens[otherGeneCol];
                    String[] otherGenes = token.split(", ");
                    for (String tmp : otherGenes)
                        genes.add(tmp);
                }
                // To keep the original sample counting
//                if (gene.equals("Unknown"))
//                    continue;
                //genes.add(gene);
            }
            fu.close();
        }
//        System.out.println("Total samples in CNVs: " + samples.size());
//        for (String sample : sampleToCNVGenes.keySet()) {
//            Set<String> genes = sampleToCNVGenes.get(sample);
//            System.out.println(sample + ": " + genes.size());
//        }
        return sampleToCNVGenes;
    }
    
    /**
     * This method is used to load glioblastoma (GBM) genes published in Science.
     * @return
     * @throws IOException
     */
    protected List<String> loadScienceGBMGenes() throws IOException {
        List<String> genes = new ArrayList<String>();
        String fileName = DIR_NAME + "GBMCanGenes.txt";
        return loadScienceGenes(genes, fileName);
    }
    
    private List<String> loadLungCancerGenes() throws IOException {
        String fileName = R3Constants.RESULT_DIR + "ScottPowers/KnownLungCancerGenes.txt";
        List<String> genes = new ArrayList<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            String tokens[] = line.split("\t");
            genes.add(tokens[0]);
        }
        fu.close();
        return genes;
    }
    
    private List<String> calculateShortestPath(HierarchicalClusterNode node1,
                                               HierarchicalClusterNode node2,
                                               Map<String, List<String>> pairToPath) {
        List<String> miniPath = null;
        for (String id1 : node1.ids) {
            for (String id2 : node2.ids) {
                String key = generateSortedKey(id1, id2);
                List<String> shortestPath = pairToPath.get(key);
                if (miniPath == null || 
                    shortestPath.size() < miniPath.size())
                    miniPath = shortestPath;
            }
        }
        return miniPath;
    }
    
    private String generateSortedKey(String id1, String id2) {
        int compare = id1.compareTo(id2);
        if (compare < 0)
            return id1 + "\t" + id2;
        return id2 + "\t" + id1;
    }
    
    /**
     * Load a map from cancer sample to mutated genes. This is based on the TCGA data format
     * table.
     * @param fileName
     * @return
     * @throws IOException
     */
    protected Map<String, Set<String>> loadSampleToMutatedGenes(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            String sample = tokens[15];
            sample = sample.substring(0, 12); // Use sample ID only
            Set<String> set = sampleToGenes.get(sample);
            if (set == null) {
                set = new HashSet<String>();
                sampleToGenes.put(sample, set);
            }
            set.add(gene);
        }
        fu.close();
        return sampleToGenes;
    }

    /**
     * Load a list of genes ranked by MSKCC method (downloaded from 
     * http://awabi.cbio.mskcc.org/tcga_gbm/jsp/index.jsp?filer_type=by_score&score_threshold=1).
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadGeneRankers() throws IOException {
        String fileName = NatureGBMAnalyzer.TCGA_GBM_DIR + "Gene_Ranker_All_Fixed.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[2];
            Double score = new Double(tokens[1]);
            geneToScore.put(gene, score);
        }
        fu.close();
        return geneToScore;
    }

    protected void hierchicallyClusterSamples(List<String> resequencedSamples,
                                              Map<String, Set<String>> sampleToAlteredGenes)
            throws IOException {
        // Clustering based on all mutations
        System.out.println("Calculate scores...");
        Map<String, Double> samplePairToScore = calculateSamplePairToDistanceScore(resequencedSamples, 
                                                                                   sampleToAlteredGenes);
//        Map<String, Double> samplePairToScore = calculateSamplePairToScoreFromShare(resequencedSamples, 
//                                                                                    sampleToAlteredGenes);
        System.out.println("hierarchical layouting...");
        Set<String> pairedSamples = InteractionUtilities.grepIDsFromInteractions(samplePairToScore.keySet());
        List<HierarchicalClusterNode> clusters = hierarchicalCluster(pairedSamples, 
                                                         samplePairToScore,
                                                         false);
        HierarchicalClusterNode root = clusters.get(0);
        outputCluster(root);
    }
    
    private Map<String, Double> calculateSamplePairToDistanceScore(List<String> resequencedSamples,
                                                                   Map<String, Set<String>> sampleToAlteredGenes) throws IOException {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> idToPartners = bfs.generateIdToPartnersMap(fis);
        // Clustering based on all mutations
        Map<String, Double> samplePairToDistance = new HashMap<String, Double>();
        Set<String> pairedSamples = new HashSet<String>();
        System.out.println("calculate shortest path...");
        for (int i = 0; i < resequencedSamples.size() - 1; i++) {
            String sample1 = resequencedSamples.get(i);
            Set<String> genes1 = sampleToAlteredGenes.get(sample1);
            genes1.retainAll(allGenes);
            if (genes1.size() == 0)
                continue;
            pairedSamples.add(sample1);
            for (int j = i + 1; j < resequencedSamples.size(); j++) {
                String sample2 = resequencedSamples.get(j);
                Set<String> genes2 = sampleToAlteredGenes.get(sample2);
                genes2.retainAll(allGenes);
                if (genes2.size() == 0)
                    continue;
                pairedSamples.add(sample2);
//                double dist = calculateAverageDistance(genes1, 
//                                                    genes2, 
//                                                    bfs, 
//                                                    idToPartners);
                double dist = calculateMinShortestPath(genes1, 
                                                       genes2, 
                                                       bfs,
                                                       idToPartners);
                samplePairToDistance.put(sample1 + "\t" + sample2,
                                         dist);
                System.out.println(sample1 + "\t" + sample2 + ": " + dist);
            }
        }
        return samplePairToDistance;
    }
    
    private double calculateMinShortestPath(Set<String> set1,
                                            Set<String> set2,
                                            BreadthFirstSearch bfs,
                                            Map<String, Set<String>> idToPartners) {
        List<String> list1 = new ArrayList<String>(set1);
        List<String> list2 = new ArrayList<String>(set2);
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
        // From list2 to list1
        for (String gene2 : list2) {
            Map<String, Integer> geneToDistances = bfs.getDistances(gene2, 
                                                                    list1,
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

    protected void filterSampleToGenes(Map<String, Set<String>> sampleToGenes, int sampleNumber) {
        CancerAnalysisUtilitites.filterSampleToGenes(sampleToGenes, sampleNumber);
    }

    protected Set<String> selectGenesInSamples(int sampleNumber, Map<String, Set<String>> sampleToGenes) {
        return CancerAnalysisUtilitites.selectGenesInSamples(sampleNumber, sampleToGenes);
    }

}
