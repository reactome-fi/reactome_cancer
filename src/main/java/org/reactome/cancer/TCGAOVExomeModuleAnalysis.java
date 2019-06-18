/*
 * Created on Jan 24, 2012
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.TransformerUtils;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.r3.TopicAnalyzer;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.JGraphTUtilities;
import org.reactome.r3.graph.JungGraphUtilities;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import edu.uci.ics.jung.algorithms.scoring.PageRankWithPriors;
import edu.uci.ics.jung.graph.Graph;

/**
 * @author gwu
 *
 */
public class TCGAOVExomeModuleAnalysis extends TCGAOvarianCancerAnalyzer {
    
    public TCGAOVExomeModuleAnalysis() {
    }
    
    @Test
    public void checkMutSigForModuleGenes() throws IOException {
        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/";
        String geneSigFile = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/" + 
                "OV.sig_genes.txt";
        String clusterFileName = dirName + "2009FISubNetwork_Sample3_Modules_092012.txt";
        String mutationFile = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/" + 
                 "OV.final_analysis_set.maf";
        int sampleCutoff = 3;
        
//        dirName = R3Constants.GBM_DIR + "FireHose/";
//        geneSigFile = dirName + "FH_GBM.Mutation_Significance.Level_4.20120725/GBM.sig_genes.txt";
//        clusterFileName = dirName + "2012FISubNetwork_Sample5_Modules_092012.txt";
//        mutationFile = dirName + "FH_GBM.Mutation_Significance.Level_4.20120725/GBM.final_analysis_set.maf";
//        sampleCutoff = 5;
        
        
        Map<String, Double> geneToPValue = loadGeneToMutSigPValue(geneSigFile);
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
        MATFileLoader matFileLoader = new MATFileLoader();
        Map<String, Set<String>> sampleToGenes = matFileLoader.loadSampleToGenes(mutationFile, false);
        filterSampleToGenes(sampleToGenes, sampleCutoff);
        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
        
        Set<String> module = clusters.get(3);
        for (String gene : module) {
            System.out.println(gene + "\t" + geneToPValue.get(gene));
        }
        
        // Want to print out all genes for Fisher test
        System.out.println("\nAll genes mutated in " + sampleCutoff);
        for (String gene : allGenes) {
            System.out.println(gene + "\t" + geneToPValue.get(gene));
        }
    }
    
    @Test
    public void generateSigGenes() throws IOException {
        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/";
        String geneSigFile = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/" + 
                             "OV.sig_genes.txt";
        Map<String, Double> geneToPvalue = loadGeneToMutSigPValue(geneSigFile);
        for (String gene : geneToPvalue.keySet()) {
            Double pvalue = geneToPvalue.get(gene);
            if (pvalue <= 0.05)
                System.out.println(gene);
        }
    }
    
    public Map<String, Double> loadGeneToMutSigPValue(String fileName) throws IOException {
        Map<String, Double> geneToPvalue = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[1];
            String pvalue = tokens[19]; // Use p-value instead of q values for our purpose.
            if (pvalue.startsWith("<"))
                pvalue = pvalue.substring(1);
            geneToPvalue.put(gene, new Double(pvalue));
        }
        fu.close();
        return geneToPvalue;
    }
    
    /**
     * Perform a network clustering using spectral partition for the TCGA OV exome
     * mutation data file using a method call.
     * @return
     * @throws IOException
     */
    private List<Set<String>> performNetworkCluster(Map<String, Set<String>> sampleToGenes,
                                                    int moduleSizeCutoff) throws IOException {
        System.out.println("Total samples: " + sampleToGenes.size());
        Set<String> genes = selectGenesInSamples(4, sampleToGenes);
        System.out.println("Total genes mutated in 4 or more samples: " + genes.size());
        SpectralPartitionNetworkCluster clustering = new SpectralPartitionNetworkCluster();
        List<Set<String>> clusters = clustering.clusterGenes(genes);
        System.out.println("Total clusters: " + clusters.size());
        int index = 0;
        for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> cluster = it.next();
            if (cluster.size() < moduleSizeCutoff)
                it.remove();
            System.out.println(index + "\t" + cluster.size());
            index ++;
        }
        return clusters;
    }
    
    /**
     * This method is used to generate a samples to network modules matrix by using
     * a weight scheme devised by Christina: each gene is assigned a weight equivalent
     * to the number of mutated samples; each module assigned a weight by adding all gene
     * weights together; a score for a sample in a module is sum of weights of genes mutated
     * in the sample divided by total the total module weight.
     * @throws Exception
     */
    @Test
    public void generateSampleToModuleWeightMatrix() throws Exception {
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        int moduleSizeCutoff = 9;
        
        List<Set<String>> clusters = performNetworkCluster(sampleToGenes, 
                                                           moduleSizeCutoff);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        // Calculate module weights
        List<Integer> clusterWeights = new ArrayList<Integer>();
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            int weight = 0;
            for (String gene : cluster) {
                weight += geneToSamples.get(gene).size();
            }
            clusterWeights.add(weight);
        }
        
        StringBuilder builder = new StringBuilder();
        generateSampleToModuleHeader(clusters, builder);
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            builder.append(sample);
            Set<String> mutatedGenes = sampleToGenes.get(sample);
            for (int i = 0; i < clusters.size(); i++) {
                builder.append("\t");
                int total = 0;
                for (String gene : clusters.get(i)) {
                    if (mutatedGenes.contains(gene)) {
                        total += geneToSamples.get(gene).size();
                    }
                }
                double score = (double) total / clusterWeights.get(i);
                builder.append(score);
            }
            System.out.println(builder.toString());
            builder.setLength(0);
        }
    }
    
    /**
     * Construct a FI sub-network, do a network clustering, and generate a binary matrix
     * for all samples.
     * @throws Exception
     */
    @Test
    public void generateSampleToModuleBinaryMatrix() throws Exception {
        int moduleSizeCutoff = 9;
        
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        List<Set<String>> clusters = performNetworkCluster(sampleToGenes,
                                                           moduleSizeCutoff);
        
        StringBuilder builder = new StringBuilder();
        generateSampleToModuleHeader(clusters, builder);
        
        // A simple binary matrix
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            builder.append(sample);
            Set<String> genes = sampleToGenes.get(sample);
            for (int i = 0; i < clusters.size(); i++) {
                boolean isShared = InteractionUtilities.isShared(genes, clusters.get(i));
                builder.append("\t");
                if (isShared)
                    builder.append("1");
                else
                    builder.append("0");
            }
            System.out.println(builder.toString());
            builder.setLength(0);
        }
    }

    private void generateSampleToModuleHeader(List<Set<String>> clusters,
                                              StringBuilder builder) {
        builder.append("Sample");
        for (int i = 0; i < clusters.size(); i++) {
            builder.append("\t").append("Module").append(i);
        }
//        System.out.println(builder.toString());
//        builder.setLength(0);
    }
    
    @Test
    public void generateSampleToMutatedGeneMatrix() throws Exception {
//        String dirName = "datasets/TCGA/OvarianCancer/FireHose/";
//        String mutationFileName = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/OV.final_analysis_set.maf"; 
//        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);

        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);

        
        filterSampleToGenes(sampleToGenes, 4);
        Set<String> totalGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
        
//        String moduleFileName = dirName + "2012FISubNetwork_Sample3_No_UBC_Modules_092012.txt";
//        List<Set<String>> modules = new NetworkClusterAnalyzer().loadNetworkClusters(moduleFileName);
//        Set<String> totalGenes = new HashSet<String>();
//        for (Set<String> module : modules)
//            totalGenes.addAll(module);
        
        System.out.println("Total samples: " + sampleToGenes.size());
        System.out.println("Total genes: " + totalGenes.size());
        
        // For PageRank scoring
        JungGraphUtilities gu = new JungGraphUtilities();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisInGenes = InteractionUtilities.getFIs(totalGenes, fis);
        // Have to work with genes in the FI network only
        totalGenes = InteractionUtilities.grepIDsFromInteractions(fisInGenes);
        Graph<String, String> graph = gu.createJungGraph(fisInGenes);
        
//        // Uniform prior
//        PageRank<String, String> ranks = new PageRank<String, String>(graph, 0.15d);
//        ranks.evaluate();
        // Prior PageRank
        Map<String, Double> prior = createPageRankPriorBasedOnSampleToGenes(sampleToGenes);
        Transformer<String, Double> transformer = TransformerUtils.mapTransformer(prior);
        PageRankWithPriors<String, String> ranks = new PageRankWithPriors<String, String>(graph, 
                                                                                          transformer,
                                                                                          0.15d);
        ranks.evaluate();
        //String fileName = dirName + "SampleToGeneMatrix_3Samples_092412.txt";
//        String fileName = dirName + "SampleToGeneMatrix_2012FISubNetwork_Sample3_No_UBC_Modules_092012.txt";
//        String fileName = dirName + "SampleToExomeGeneBinary_2012_4Sample_071713.txt";
//        String fileName = dirName + "SampleToExomeGenePageRank_2012_4Sample_071713.txt";
        String fileName = dirName + "SampleToExomeGenePageRankPrior_2012_4Sample_071713.txt";
        
        StringBuilder builder = new StringBuilder();
        fu.setOutput(fileName);
        builder.append("SampleId");
        
        List<String> geneList = new ArrayList<String>(totalGenes);
        Collections.sort(geneList);
        for (String gene : geneList) {
            builder.append("\t").append(gene);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            builder.append(sample);
            for (String gene : geneList) {
                //builder.append("\t").append(genes.contains(gene) ? "1" : "0");
                builder.append("\t");
                if (genes.contains(gene))
                    builder.append(ranks.getVertexScore(gene));
                else
                    builder.append("0");
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkOverlappingBetweenOldModulesAndNewModule() throws Exception {
        String moduleGenesFile = OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> oldModules = fu.loadInteractions(moduleGenesFile);
        System.out.println("Old modules: " + oldModules.size());
//        String moduleFile = OVARIAN_DIR_NAME + "data_published_in_nature/OVNetowrkModules2012_091812.txt";
        String moduleFile = OVARIAN_DIR_NAME + "FireHose/2009FISubNetwork_Sample3_Modules_092012.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> modules = clusterAnalyzer.loadNetworkClusters(moduleFile);
        Set<String> newModule = modules.get(3);
        System.out.println("New module: " + newModule.size());
        Set<String> shared = InteractionUtilities.getShared(oldModules, newModule);
        System.out.println("Shared genes between 2009 modules and 2012 module: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_PROTEIN_IDS, 
                                                                    oldModules.size(),
                                                                    newModule.size(),
                                                                    shared.size());
        System.out.println("P-value: " + pvalue);
    }
    
    @Test
    public void checkOverlappingBetweenSuperpcAndModulesForPathways() throws Exception {
        Set<String> superpcGenes = new OvarianMCLGeneExpModuleAnalyzer().getSelectedGenesFromSuperpc();
        System.out.println("Superpc Genes: " + superpcGenes.size());
        String moduleGenesFile = OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> moduleGenes = fu.loadInteractions(moduleGenesFile);
        TopicAnalyzer pathwayAnalyzer = new TopicAnalyzer();
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.05d);
        System.out.println("Annotate superpc genes:");
        annotator.annotateGenesWithFDR(superpcGenes, AnnotationType.Pathway);
        System.out.println("\nAnnotate module genes:");
        annotator.annotateGenesWithFDR(moduleGenes, AnnotationType.Pathway);
    }
    
    /**
     * Check a list of random gene sets that are strong correlated with TCGA
     * OV patient overall survivals.
     * @throws Exception
     */
    @Test
    public void checkBestRandomGeneSets() throws Exception {
        Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        filterSampleToGenes(sampleToGenes, 3);
        Set<String> geneSet = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Total genes: " + geneSet.size());

        Set<String> allFIs = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(allFIs);
        Set<String> fis = InteractionUtilities.getFIs(geneSet, allFIs);
        
        List<Set<String>> modules = new SpectralPartitionNetworkCluster().cluster(fis);
        System.out.println("Total modules: " + modules.size());
        
        Set<String> pickedModule = modules.get(3);
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.25d);
        System.out.println("Picked TCGA OV Exome Module:");
        annotator.annotateGenesWithFDR(pickedModule, AnnotationType.Pathway);
        System.out.println();
        
        String fileName = OVARIAN_DIR_NAME + "FireHose/BestRandomGeneSets_051413.txt";
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        List<List<String>> randomGeneSets = new ArrayList<List<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(",");
            List<String> genes = Arrays.asList(tokens);
            System.out.println("Random module " + index);
            index ++;
            annotator.annotateGenesWithFDR(genes, AnnotationType.Pathway);
            System.out.println();
            randomGeneSets.add(genes);
        }
        fu.close();
        
        System.out.println("\nShortest path");
        // Check shortest path
        BreadthFirstSearch searcher = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = searcher.generateIdToPartnersMap(fis);
        NatureGBMAnalyzer analyzer = new NatureGBMAnalyzer();
        for (int i = 0; i < randomGeneSets.size(); i++) {
            double distance = analyzer.calculateAverageDistance(pickedModule,
                                                                new HashSet<String>(randomGeneSets.get(i)),
                                                                searcher, 
                                                                idToPartners);
            System.out.println("Random gene set " + i + ": " + distance);
        }

    }
    
    @Test
    public void checkOverlappingBetweenSuperpcAndModules() throws Exception {
        Set<String> superpcGenes = new OvarianMCLGeneExpModuleAnalyzer().getSelectedGenesFromSuperpc();
        System.out.println("Superpc Genes: " + superpcGenes.size());
        String moduleGenesFile = OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> moduleGenes = fu.loadInteractions(moduleGenesFile);
        Set<String> shared = InteractionUtilities.getShared(superpcGenes, moduleGenes);
        System.out.println("Module 6 + 7 genes: " + moduleGenes.size());
        System.out.println("Shared: " + shared.size());
        // Check how many FIs between these two gene sets
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        shared = InteractionUtilities.getShared(allGenes, superpcGenes);
        System.out.println("Superpc Genes in Biggest component: " + shared.size());
        shared = InteractionUtilities.getShared(allGenes, moduleGenes);
        System.out.println("Module genes in biggest component: " + shared.size());
        // Make sure all genes in the connected component
        int numberOfFIsBetween = countFIsNumberBetween(superpcGenes,
                                                       moduleGenes, 
                                                       fis);
        System.out.println("Total FIs: " + numberOfFIsBetween);
        // Check shortest path
        BreadthFirstSearch searcher = new BreadthFirstSearch();
        Map<String, Set<String>> idToPartners = searcher.generateIdToPartnersMap(fis);
        NatureGBMAnalyzer analyzer = new NatureGBMAnalyzer();
        double distance = analyzer.calculateAverageDistance(superpcGenes,
                                                            moduleGenes,
                                                            searcher, 
                                                            idToPartners);
        System.out.println("Average distance: " + distance);
        // Just a simple permutation test
        int permutate = 1000;
        int size = superpcGenes.size() + moduleGenes.size();
        RandomData randomizer = new RandomDataImpl();
        List<Integer> counts = new ArrayList<Integer>();
        List<Double> randomAverageDists = new ArrayList<Double>();
        for (int i = 0; i < permutate; i++) {
            System.out.println("Permutate " + i + "...");
            Set<String> randomSet = org.reactome.r3.util.MathUtilities.randomSampling(allGenes,
                                                                                      size,
                                                                                      randomizer);
            // Need to split it into two sets
            Set<String> firstSet = MathUtilities.randomSampling(randomSet,
                                                                superpcGenes.size(),
                                                                randomizer);
            // Used as the second set
            randomSet.removeAll(firstSet);
            int randomFIsBetween = countFIsNumberBetween(firstSet, randomSet, fis);
//            System.out.println(randomFIsBetween);
            counts.add(randomFIsBetween);
            double randomDist = analyzer.calculateAverageDistance(firstSet, randomSet, searcher, idToPartners);
            randomAverageDists.add(randomDist);
        }
//        Collections.sort(counts);
        System.out.println("\nCounts from random tests:");
        int targetIndex = 0;
        for (Integer count : counts) {
//            System.out.println(count);
            if (count <= numberOfFIsBetween)
                targetIndex ++;
        }
        System.out.println("Target index: " + targetIndex);
        System.out.println("pvalue (for FIs): " + targetIndex / (double) permutate);
//        Collections.sort(randomAverageDists);
        // Want to output
        String outFileName = OVARIAN_DIR_NAME + "data_published_in_nature/AverageShortestPathRandom1000_041912.txt";
        fu.setOutput(outFileName);
        for (int i = 0; i < randomAverageDists.size(); i++) {
//        for (Double dist : randomAverageDists)
            Double dist = randomAverageDists.get(i);
            Integer count = counts.get(i);
            fu.printLine(dist + "\t" + count);
        }
        fu.close();
    }

    private int countFIsNumberBetween(Set<String> superpcGenes,
                                      Set<String> moduleGenes, 
                                      Set<String> fis) {
        int numberOfFIsBetween = 0;
        for (String gene1 : superpcGenes) {
            for (String gene2 : moduleGenes) {
                if (fis.contains(gene1 + "\t" + gene2) ||
                    fis.contains(gene2 + "\t" + gene1))
                    numberOfFIsBetween ++;
            }
        }
        return numberOfFIsBetween;
    }
    
    @Test
    public void doMutalExclusityAnalysis() throws Exception {
        Map<String, Set<String>> sampleToMutatedGenes = loadSampleToNonSynonymousMutatedGenes(false);
        System.out.println("Total samples: " + sampleToMutatedGenes.size());
        Set<String> totalMutatedGenes = InteractionUtilities.grepAllGenes(sampleToMutatedGenes);
        System.out.println("Total mutated genes: " + totalMutatedGenes.size());
//        filterSampleToGenes(sampleToMutatedGenes, 3);
//        totalMutatedGenes = InteractionUtilities.grepAllGenes(sampleToMutatedGenes);
//        System.out.println("Total mutated genes (in 3 or more samples): " + totalMutatedGenes.size());
//        // Load FIs
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
//        Set<String> allFIGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        totalMutatedGenes.retainAll(allFIGenes);
//        System.out.println("In FIs: " + totalMutatedGenes.size());
//        // Make sure all genes linked
//        Set<String> fisInMutatedGenes = InteractionUtilities.getFIs(totalMutatedGenes, fis);
//        totalMutatedGenes = InteractionUtilities.grepIDsFromInteractions(fisInMutatedGenes);
//        System.out.println("Linked: " + totalMutatedGenes.size());
        
        String moduleGenesFile = OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> moduleGenes = fu.loadInteractions(moduleGenesFile);
        System.out.println("\nTotal module genes: " + moduleGenes.size() + ", " + moduleGenes);
//        String text = "KCNA3, DYNC1I1, PCNT, GRIA4, CACNA1F";
//        String[] tokens = text.split(", ");
//        Set<String> toBeRemoved = new HashSet<String>();
//        for (String token : tokens)
//            toBeRemoved.add(token);
//        moduleGenes.removeAll(toBeRemoved);
        // Generate Sample To mutated genes
        FisherExact fisherTest = new FisherExact(totalMutatedGenes.size());
        System.out.println("Sample\tTotalMutated\tNumberOfGenes\tMutatedGenes\tP_value");
        for (String sample : sampleToMutatedGenes.keySet()) {
            Set<String> mutatedGenes = sampleToMutatedGenes.get(sample);
            mutatedGenes.retainAll(totalMutatedGenes);
            Set<String> shared = InteractionUtilities.getShared(mutatedGenes, moduleGenes);
            if (shared.size() == 0)
                continue;
            List<String> geneList = new ArrayList<String>(shared);
            Collections.sort(geneList);
            String temp = geneList.toString();
//            double pvalue = MathUtilities.calculateHypergeometricPValue(totalMutatedGenes.size(), 
//                                                                        mutatedGenes.size(), 
//                                                                        moduleGenes.size(),
//                                                                        shared.size());
            double pvalue = fisherTest.getRightTailedP(shared.size(),
                                                       moduleGenes.size() - shared.size(), 
                                                       mutatedGenes.size() - shared.size(), 
                                                       totalMutatedGenes.size() - moduleGenes.size() - mutatedGenes.size() + shared.size());
            System.out.println(sample + "\t" + 
                               mutatedGenes.size() + "\t" + 
                               shared.size() + "\t" + 
                               temp.substring(1, temp.length() - 1) + "\t" +
                               String.format("%.4f",pvalue));
        }
        
        // Generate a binary matrix
        System.out.println("\nBinary matrix");
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        List<String> moduleGeneList = new ArrayList<String>(moduleGenes);
        Collections.sort(moduleGeneList);
        for (String gene : moduleGeneList)
            builder.append("\t").append(gene);
        System.out.println(builder.toString());
        builder.setLength(0);
        for (String sample : sampleToMutatedGenes.keySet()) {
            builder.append(sample);
            Set<String> mutatedGenes = sampleToMutatedGenes.get(sample);
            for (String gene : moduleGeneList) {
                builder.append("\t").append(mutatedGenes.contains(gene) ? 1 : 0);
            }
            System.out.println(builder.toString());
            builder.setLength(0);
        }
        
//        // Generate a binary graph
//        System.out.println("\nBinary graph:");
//        for (String sample : sampleToMutatedGenes.keySet()) {
//            Set<String> mutatedGeSet = sampleToMutatedGenes.get(sample);
//            mutatedGeSet = new HashSet<String>(mutatedGeSet);
//            mutatedGeSet.retainAll(moduleGenes);
//            for (String gene : mutatedGeSet)
//                System.out.println(sample + "\tSampleToMutation\t" + gene);
//        }
//        Set<String> fisInModule = InteractionUtilities.getFIs(moduleGenes, fis);
//        for (String fi : fisInModule) {
//            String[] tokens = fi.split("\t");
//            System.out.println(tokens[0] + "\tFI\t" + tokens[1]);
//        }
//        System.out.println("ID\tType");
//        for (String gene : moduleGenes)
//            System.out.println(gene + "\tGene");
    }
    
    /**
     * This method is used to prune the module by using an augumented FI network by attaching
     * samples to the modules.
     * @throws Exception
     */
    @Test
    public void pruneModule() throws Exception {
        // Load FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> allFIGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        
        String moduleGenesFile = OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> moduleGenes = fu.loadInteractions(moduleGenesFile);
        System.out.println("\nTotal module genes: " + moduleGenes.size() + ", " + moduleGenes);
        
        // Construct a module graph
        Set<String> fisInModule = InteractionUtilities.getFIs(moduleGenes, fis);
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        SimpleGraph<String, DefaultEdge> moduleGraph = (SimpleGraph<String, DefaultEdge>) JGraphTUtilities.createGraph(moduleGenes, fisInModule);
        // Construct an augmented graph
        Map<String, Set<String>> sampleToMutatedGenes = loadSampleToNonSynonymousMutatedGenes(false);
        Set<String> selectedSamples = new HashSet<String>();
        Set<String> sampleToGene = new HashSet<String>();
        for (String sample : sampleToMutatedGenes.keySet()) {
            Set<String> mutatedGenes = sampleToMutatedGenes.get(sample);
            Set<String> sharedGenes = InteractionUtilities.getShared(mutatedGenes, moduleGenes);
            if (sharedGenes.size() == 0)
                continue;
            for (String gene : sharedGenes) {
                sampleToGene.add(sample + "\t" + gene);
            }
            selectedSamples.add(sample);
        }
        Set<String> augumentNodes = new HashSet<String>(moduleGenes);
        augumentNodes.addAll(selectedSamples);
        Set<String> augumentEdges = new HashSet<String>(fisInModule);
        augumentEdges.addAll(sampleToGene);
        SimpleGraph<String, DefaultEdge> augumentedGraph = (SimpleGraph<String, DefaultEdge>)JGraphTUtilities.createGraph(augumentNodes, augumentEdges);
        
        // Start pruning
        for (String gene : moduleGenes) {
            boolean canPrune = canPrune(moduleGraph, augumentedGraph, gene);
            if (canPrune) {
                System.out.println(gene);
                pruneGene(moduleGraph, augumentedGraph, gene, moduleGenes, "\t");
//                break;
            }
        }
    }
    
    private void pruneGene(SimpleGraph<String, DefaultEdge> moduleGraph,
                           SimpleGraph<String, DefaultEdge> auguGraph,
                           String prunedGene,
                           Set<String> checkingGenes,
                           String prefix) {
        checkingGenes = new HashSet<String>(checkingGenes);
        checkingGenes.remove(prunedGene);
        moduleGraph = (SimpleGraph<String, DefaultEdge>) moduleGraph.clone();
        moduleGraph.removeVertex(prunedGene);
        auguGraph = (SimpleGraph<String, DefaultEdge>) auguGraph.clone();
        auguGraph.removeVertex(prunedGene);
        for (String gene : checkingGenes) {
            boolean canPrune = canPrune(moduleGraph, auguGraph, gene);
            if (canPrune) {
                System.out.println(prefix + "" + gene);
                pruneGene(moduleGraph, auguGraph, gene, checkingGenes, prefix + "\t");
//                break;
            }
        }
    }
    
    private boolean canPrune(SimpleGraph<String, DefaultEdge> moduleGraph,
                             SimpleGraph<String, DefaultEdge> auguGraph,
                             String gene) {
        SimpleGraph<String, DefaultEdge> clone = (SimpleGraph<String, DefaultEdge>) moduleGraph.clone();
        clone.removeVertex(gene);
        int compNumber = getComponentNumber(clone);
        if (compNumber > 1)
            return false;
        clone = (SimpleGraph<String, DefaultEdge>) auguGraph.clone();
        clone.removeVertex(gene);
        compNumber = getComponentNumber(clone);
        if (compNumber == 1)
            return true;
        return false;
    }
    
    private int getComponentNumber(SimpleGraph<String, DefaultEdge> graph) {
        ConnectivityInspector<String, DefaultEdge> inspector = new ConnectivityInspector<String, DefaultEdge>(graph);
        List<Set<String>> components = inspector.connectedSets();
        return components.size();
    }
    
    private Map<String, Set<String>> getOVPathwayToGenes() throws Exception {
        // Load pathway genes
        AnnotationHelper annotationHelper = new AnnotationHelper();
        annotationHelper.setProteinNameToPathwayFile(R3Constants.RESULT_DIR + "ProteinNameToTopics072512.txt");
        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToPathwaysMap();
        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
        //        for (String pathway : pathwayToGenes.keySet()) {
        //            Set<String> genes = pathwayToGenes.get(pathway);
        //            System.out.println(pathway + ": " + genes.size());
        //        }
        //      String[] pathways = new String[] {
        //              "Calcium signaling pathway(K)",
        //              "Arrhythmogenic right ventricular cardiomyopathy (ARVC)(K)",
        //              "Cardiac muscle contraction(K)",
        //              "Hypertrophic cardiomyopathy (HCM)(K)",
        //              "MAPK signaling pathway(K)",
        //              "Dilated cardiomyopathy(K)",
        //              "Retrograde endocannabinoid signaling(K)",
        //              "Heterotrimeric G-protein signaling pathway-Gq alpha and Go alpha mediated pathway(P)",
        //              "Serotonergic synapse(K)",
        //              "5HT2 type receptor mediated signaling pathway(P)",
        //              "Beta1 adrenergic receptor signaling pathway(P)",
        //              "GABAergic synapse(K)",
        //              "Cholinergic synapse(K)",
        //              "Type II diabetes mellitus(K)",
        //              "Glutamatergic synapse(K)",
        //              "Alzheimer's disease(K)",
        //              "GnRH signaling pathway(K)",
        //              "Vascular smooth muscle contraction(K)",
        //              "Integration of energy metabolism(R)",
        //              "Heterotrimeric G-protein signaling pathway-Gi alpha and Gs alpha mediated pathway(P)",
        //              "Long-term depression(K)",
        //              "Transmission across Chemical Synapses(R)"
        //      };
        //      for (String pathway : pathways) {
        //          Set<String> genes = pathwayToGenes.get(pathway);
        //          System.out.println("\n" + pathway + ": " + genes.size() + " genes");
        //          // Load mutated genes
        //          Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        //          String survivalScript = "/Users/gwu/Documents/EclipseWorkspace/caBigR3WebApp/WebContent/WEB-INF/CGISurvivalAnalysis.R";
        //          SurvivalAnalysisHelper survivalHelper = new SurvivalAnalysisHelper();
        //          survivalHelper.setrScript(survivalScript);
        //          survivalHelper.setTempDirName("tmp");
        //          String clinFileName = OVARIAN_DIR_NAME + "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt";
        //          
        //          doSurvivalAnalysisForMutationModule(survivalHelper, 
        //                                              sampleToGenes,
        //                                              genes, 
        //                                              clinFileName,
        //                                              pathway + "_");
        //      }
        //      
        String[] pathways = new String[] {
                "Calcium signaling pathway(K)",
                //              "Arrhythmogenic right ventricular cardiomyopathy (ARVC)(K)",
                //              "Cardiac muscle contraction(K)",
                //              "Hypertrophic cardiomyopathy (HCM)(K)",
                //              "Dilated cardiomyopathy(K)",
        };
        Map<String, Set<String>> rtn = new HashMap<String, Set<String>>();
        for (String pathway : pathways) {
            Set<String> genes = pathwayToGenes.get(pathway);
            rtn.put(pathway, genes);
        }
        return rtn;
    }
    
    private Map<String, Set<String>> getPathwayToGenes() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost", 
                                            "gk_central_102312",
                                            "root", 
                                            "macmysql01");
        Map<String, Set<String>> pathwayToGenes = new HashMap<String, Set<String>>();
        // Get Reactome Axon Guidance and its sub-pathway
        Long dbId = 422475L;
//        // Platelet activation, signaling and aggregation
//        dbId = 76002L; 
        GKInstance pathway = dba.fetchInstance(dbId);
        Set<String> genes = getGenes(pathway);
        pathwayToGenes.put(pathway.getDisplayName(), genes);
        List<GKInstance> subPathways = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        for (GKInstance subPathway : subPathways) {
            genes = getGenes(subPathway);
            pathwayToGenes.put(subPathway.getDisplayName(), genes);
        }
        return pathwayToGenes;
    }
    
    private Set<String> getGenes(GKInstance pathway) throws Exception {
        Set<GKInstance> refEntities = InstanceUtilities.grepRefPepSeqsFromPathway(pathway);
        Set<String> genes = new HashSet<String>();
        for (GKInstance refEnt : refEntities) {
            if (refEnt.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)) {
                String geneName = (String) refEnt.getAttributeValue(ReactomeJavaConstants.geneName);
                if (geneName != null)
                    genes.add(geneName);
            }
        }
        return genes;
    }
    
    /**
     * Do surival analysis for a list of modules stored in a file.
     * @throws Exception
     */
    @Test
    public void checkSurvialForModuleFile() throws Exception {
        List<String> cancers = new ArrayList<String>();
        List<Map<String, Set<String>>> sampleToGeneList = new ArrayList<Map<String,Set<String>>>();
        List<String> clinFileNames = new ArrayList<String>();
        loadTCGAFHDataSets(cancers, sampleToGeneList, clinFileNames);
        
        String moduleFile = R3Constants.OVARIAN_DIR_NAME + "FireHose/Clusters2009OfGenesWithMutSigPValue_0.05_121212.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> modules = clusterAnalyzer.loadNetworkClusters(moduleFile);
        
        // Just a simple test
////        String[] genes = new String[]{"ELAVL2", "ELAVL3", "HSF4", "RBM39", "REST"};
//        String[] genes = new String[]{"DNAH3", "DNAI2", "DNAH7", "DNAH8", "DNAH5"};
//        modules.clear();
//        modules.add(new HashSet<String>(Arrays.asList(genes)));
        
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        for (int i = 0; i < modules.size(); i++) {
            Set<String> module = modules.get(i);
            if (module.size() == 2)
                continue;
            System.out.println("Module " + i);
            doSurvivalAnalysisForMutationModule(survivalHelper, 
                                                sampleToGeneList.get(0),
                                                module, 
                                                clinFileNames.get(0), 
                                                "Module_" + i + "_");
        }
    }
    
    /**
     * Test mutations in a specific pathway.
     * @throws Exception
     */
    @Test
    public void testMutationsInPathway() throws Exception {
        // Load pathway genes
//        AnnotationHelper annotationHelper = new AnnotationHelper();
//        annotationHelper.setProteinNameToPathwayFile(R3Constants.RESULT_DIR + "ProteinNameToTopics072512.txt");
//        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToPathwaysMap();
//        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
//        //        for (String pathway : pathwayToGenes.keySet()) {
//        //            Set<String> genes = pathwayToGenes.get(pathway);
//        //            System.out.println(pathway + ": " + genes.size());
//        //        }
////        String[] pathways = new String[] {
////                "Calcium signaling pathway(K)",
////                "Arrhythmogenic right ventricular cardiomyopathy (ARVC)(K)",
////                "Cardiac muscle contraction(K)",
////                "Hypertrophic cardiomyopathy (HCM)(K)",
////                "MAPK signaling pathway(K)",
////                "Dilated cardiomyopathy(K)",
////                "Retrograde endocannabinoid signaling(K)",
////                "Heterotrimeric G-protein signaling pathway-Gq alpha and Go alpha mediated pathway(P)",
////                "Serotonergic synapse(K)",
////                "5HT2 type receptor mediated signaling pathway(P)",
////                "Beta1 adrenergic receptor signaling pathway(P)",
////                "GABAergic synapse(K)",
////                "Cholinergic synapse(K)",
////                "Type II diabetes mellitus(K)",
////                "Glutamatergic synapse(K)",
////                "Alzheimer's disease(K)",
////                "GnRH signaling pathway(K)",
////                "Vascular smooth muscle contraction(K)",
////                "Integration of energy metabolism(R)",
////                "Heterotrimeric G-protein signaling pathway-Gi alpha and Gs alpha mediated pathway(P)",
////                "Long-term depression(K)",
////                "Transmission across Chemical Synapses(R)"
////        };
////        for (String pathway : pathways) {
////            Set<String> genes = pathwayToGenes.get(pathway);
////            System.out.println("\n" + pathway + ": " + genes.size() + " genes");
////            // Load mutated genes
////            Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
////            String survivalScript = "/Users/gwu/Documents/EclipseWorkspace/caBigR3WebApp/WebContent/WEB-INF/CGISurvivalAnalysis.R";
////            SurvivalAnalysisHelper survivalHelper = new SurvivalAnalysisHelper();
////            survivalHelper.setrScript(survivalScript);
////            survivalHelper.setTempDirName("tmp");
////            String clinFileName = OVARIAN_DIR_NAME + "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt";
////            
////            doSurvivalAnalysisForMutationModule(survivalHelper, 
////                                                sampleToGenes,
////                                                genes, 
////                                                clinFileName,
////                                                pathway + "_");
////        }
////        
//        String[] pathways = new String[] {
////                "Calcium signaling pathway(K)",
////                "Arrhythmogenic right ventricular cardiomyopathy (ARVC)(K)",
////                "Cardiac muscle contraction(K)",
////                "Hypertrophic cardiomyopathy (HCM)(K)",
////                "Dilated cardiomyopathy(K)",
//                "Axon guidance(R)"
//        };
//        Set<String> shared = new HashSet<String>();
//        for (String pathway : pathways) {
//            Set<String> genes = pathwayToGenes.get(pathway);
//            shared.addAll(genes);
////            if (shared.size() == 0)
////                shared = new HashSet<String>(genes);
////            else
////                shared = InteractionUtilities.getShared(shared, genes);
//        }
//        System.out.println("Total genes in calcilum pathway: " + shared.size());
//////        for (String gene : shared)
////            System.out.println(gene);
        
        List<String> cancers = new ArrayList<String>();        
        List<Map<String, Set<String>>> sampleToGenesList = new ArrayList<Map<String,Set<String>>>();
        List<String> clinFileNames = new ArrayList<String>();
        loadTCGAFHDataSets(cancers, sampleToGenesList, clinFileNames);
        
//        
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        
//        System.out.println("\nUsing the whole KEGG Calcium Signaling Pathway:");
//        Map<String, Set<String>> pathwayToGenes = getPathwayToGenes();
        
        Map<String, Set<String>> pathwayToGenes = getOVPathwayToGenes();
        
        for (String pathway : pathwayToGenes.keySet()) {
            Set<String> shared = pathwayToGenes.get(pathway);
            System.out.println("\nPathway \"" + pathway + "\" with " + shared.size() + " genes:");
            for (int i = 0; i < cancers.size(); i++) {
                String cancer = cancers.get(i);
                System.out.println("\nAnalyzing TCGA " + cancer + "...");
                Map<String, Set<String>> sampleToGenes = sampleToGenesList.get(i);
                String clinFileName = clinFileNames.get(i);
                doSurvivalAnalysisForMutationModule(survivalHelper, 
                                                    sampleToGenes,
                                                    shared, 
                                                    clinFileName,
                                                    pathway + "_" + cancer);
            }
        }
//        
        String moduleFile = OVARIAN_DIR_NAME + "data_published_in_nature/OVNetowrkModules2012_091812.txt";
        moduleFile = OVARIAN_DIR_NAME + "FireHose/2009FISubNetwork_Sample3_Modules_092012.txt";
//        moduleFile = OVARIAN_DIR_NAME + "NetworkModulesCalciumSignalingPathway.txt";
//        moduleFile = "/Users/gwu/Desktop/CalciumSignalingPathwayClusters.txt";
//        moduleFile = OVARIAN_DIR_NAME + "data_published_in_nature/NetworkModules_101311.txt";
//        moduleFile = R3Constants.GBM_DIR + "Firehose/2012FISubNetwork_Sample5_Modules_092012.txt";
        int neededModule = 0;
        neededModule = 3;
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> modules = clusterAnalyzer.loadNetworkClusters(moduleFile);
//        Set<String> moduleGenes = new HashSet<String>();
//        moduleGenes.addAll(modules.get(6));
//        moduleGenes.addAll(modules.get(7));
//        System.out.println("Total module genes: " + moduleGenes.size());
//      
//        System.out.println("\nUsing the first network module in the KEGG Calcium Signaling Pathway:");
        System.out.println("\nUsing the Exome FI network module:");
        for (int i = 0; i < cancers.size(); i++) {
            String cancer = cancers.get(i);
            System.out.println("\nAnalyze TCGA " + cancer + "...");
            Map<String, Set<String>> sampleToGenes = sampleToGenesList.get(i);
            String clinFileName = clinFileNames.get(i);
            for (int j = 0; j < modules.size(); j++) {
                if (j != neededModule)
                    continue;
                System.out.println("\nModule " + j);
                Set<String> module = modules.get(j);
                doSurvivalAnalysisForMutationModule(survivalHelper, 
                                                    sampleToGenes,
                                                    module, 
                                                    clinFileName,
                                                    cancer + "_Module_" + j);
//                break;
            }
        }
    }
    
    private Map<String, Double> createPageRankPriorBasedOnSampleToGenes(Map<String, Set<String>> sampleToGenes) {
        Map<String, Double> geneToPrior = new HashMap<String, Double>();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Map<String, Double> geneToRatio = new HashMap<String, Double>();
        int totalSample = sampleToGenes.size();
        double total = 0.0d;
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            Double ratio = (double) (samples.size()) / totalSample;
            geneToRatio.put(gene, ratio);
            total += ratio;
        }
        for (String gene : geneToSamples.keySet()) {
            Double ratio = geneToRatio.get(gene);
            double prior = ratio / total;
            geneToPrior.put(gene, prior);
        }
//        for (String gene : geneToPrior.keySet())
//            System.out.println(gene + "\t" + geneToPrior.get(gene));
        return geneToPrior;
    }
    
    /**
     * This method is used to generate a matrix from sample to network module to gene degrees
     * in the modules. The score calculation is similar to gene connection degree based method except
     * page rank is used in this method.
     */
    @Test
    public void generateSampleToModuleGenePageRankMatrix() throws Exception {
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        
        int moduleSizeCutoff = 9;
        String output = dirName + "SampleToExomeModuleGenePageRankPrior_2012_4Sample_9Genes_071713.txt";
        
        List<Set<String>> clusters = performNetworkCluster(sampleToGenes, 
                                                           moduleSizeCutoff);
        // The following code is used for using prior pageranking algorithm
        // Create priors based on mutated samples
        Map<String, Double> geneToSamplePrior = createPageRankPriorBasedOnSampleToGenes(sampleToGenes);
        Transformer<String, Double> geneToSamplePriorTransformer = TransformerUtils.mapTransformer(geneToSamplePrior);
        
        // Get module weights for each cluster
        Map<String, Double> geneToPageRank = new HashMap<String, Double>();
        List<Double> clusterWeights = new ArrayList<Double>();
        // Need to load all interactions first
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        JungGraphUtilities gu = new JungGraphUtilities();
        double alpha = 0.15d;
        
        for (Set<String> cluster : clusters) {
            Set<String> clusterInterctions = InteractionUtilities.getFIs(cluster, fis);
            Graph<String, String> graph = gu.createJungGraph(clusterInterctions);
            
            // Uniform prior
//            PageRank<String, String> ranks = new PageRank<String, String>(graph, alpha);
//            ranks.evaluate();
            // Prior page rank
            PageRankWithPriors<String, String> ranks = new PageRankWithPriors<String, String>(graph, 
                    geneToSamplePriorTransformer, 
                    alpha);
            ranks.evaluate();
            
            Double total = 0.0d;
            for (String gene : cluster) {
                double pageRank = ranks.getVertexScore(gene);
                total += pageRank;
                geneToPageRank.put(gene, pageRank);
            }
            clusterWeights.add(total);
        }
        
        StringBuilder builder = new StringBuilder();
        generateSampleToModuleHeader(clusters, builder);
        
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            builder.append("\n");
            builder.append(sample);
            Set<String> mutatedGenes = sampleToGenes.get(sample);
            int index = 0;
            for (Set<String> cluster : clusters) {
                Set<String> shared = InteractionUtilities.getShared(cluster, mutatedGenes);
                double total = 0.0d;
                for (String gene : shared) {
                    total += geneToPageRank.get(gene);
                }
                builder.append("\t").append(total / clusterWeights.get(index));
                index ++;
            }
//            System.out.println(builder.toString());
//            builder.setLength(0);
        }
        fu.setOutput(output);
        fu.printLine(builder.toString());
        fu.close();
    }
    
    /**
     * This method is used to generate a matrix from sample to network module to gene degrees
     * in the modules. The score calculation is similar to gene connection degree based method except
     * centrality is used in this method.
     */
    @Test
    public void generateSampleToModuleGeneCentralityMatrix() throws Exception {
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        int moduleSizeCutoff = 9;
        String output = dirName + "SampleToExomeModuleGeneCentrality_2012_4Sample_9Genes_071713.txt";
        
        List<Set<String>> clusters = performNetworkCluster(sampleToGenes, 
                                                           moduleSizeCutoff);
        // Get module weights for each cluster
        Map<String, Double> geneToCentrality = new HashMap<String, Double>();
        List<Double> clusterWeights = new ArrayList<Double>();
        // Need to load all interactions first
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        JungGraphUtilities gu = new JungGraphUtilities();
        for (Set<String> cluster : clusters) {
            Set<String> clusterInterctions = InteractionUtilities.getFIs(cluster, fis);
            Map<String, Double> clusterGeneToCentrality = gu.calculateCentralities(cluster, clusterInterctions);
            geneToCentrality.putAll(clusterGeneToCentrality);
            Double total = 0.0d;
            for (String gene : clusterGeneToCentrality.keySet())
                total += clusterGeneToCentrality.get(gene);
//            if (total.equals(0.0d))
//                System.out.println("Total is 0.0: " + cluster);
            clusterWeights.add(total);
        }
        
        StringBuilder builder = new StringBuilder();
        generateSampleToModuleHeader(clusters, builder);
        
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            builder.append("\n");
            builder.append(sample);
            Set<String> mutatedGenes = sampleToGenes.get(sample);
            int index = 0;
            for (Set<String> cluster : clusters) {
                Set<String> shared = InteractionUtilities.getShared(cluster, mutatedGenes);
                double total = 0.0d;
                for (String gene : shared) {
                    total += geneToCentrality.get(gene);
                }
                builder.append("\t").append(total / clusterWeights.get(index));
                index ++;
            }
//            System.out.println(builder.toString());
//            builder.setLength(0);
        }
        fu.setOutput(output);
        fu.printLine(builder.toString());
        fu.close();
    }
    
    /**
     * This method is used to generate a matrix from sample to network module to gene degrees
     * in the modules. The score for degree is calculated as followings:
     * 1). Each module is treated as a graph. Each vertex in the graph should have a degree.
     * 2). The module graph has a total degree which is a sum of all vertex degrees.
     * 3). The score for a sample in a module is calculated as: sum_of_degrees_for_genes_mutated_in_Sample/module_degree
     * @throws Exception
     */
    @Test
    public void generateSampleToModuleGeneDegreeMatrix() throws Exception {
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/052013/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf"; 
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        int moduleSizeCutoff = 3;
        String output = dirName + "SampleToExomeModuleGeneDegree_2012_4Sample_3Genes_071713.txt";
        
        List<Set<String>> clusters = performNetworkCluster(sampleToGenes, 
                                                           moduleSizeCutoff);
        // Get module weights for each cluster
        Map<String, Integer> geneToDegree = new HashMap<String, Integer>();
        List<Integer> clusterWeights = new ArrayList<Integer>();
        // Need to load all interactions first
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        for (Set<String> cluster : clusters) {
            Set<String> clusterInterctions = InteractionUtilities.getFIs(cluster, fis);
            Map<String, Integer> clusterGeneToDegree = InteractionUtilities.generateProteinToDegree(clusterInterctions);
            geneToDegree.putAll(clusterGeneToDegree);
            int total = 0;
            for (String gene : clusterGeneToDegree.keySet())
                total += clusterGeneToDegree.get(gene);
            clusterWeights.add(total);
        }
        
        StringBuilder builder = new StringBuilder();
        generateSampleToModuleHeader(clusters, builder);
        
        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            builder.append("\n");
            builder.append(sample);
            Set<String> mutatedGenes = sampleToGenes.get(sample);
            int index = 0;
            for (Set<String> cluster : clusters) {
                Set<String> shared = InteractionUtilities.getShared(cluster, mutatedGenes);
                double total = 0.0d;
                for (String gene : shared) {
                    total += geneToDegree.get(gene);
                }
                builder.append("\t").append(total / clusterWeights.get(index));
                index ++;
            }
//            System.out.println(builder.toString());
//            builder.setLength(0);
        }
        fu.setOutput(output);
        fu.printLine(builder.toString());
        fu.close();
    }

    /**
     * This method is used to generate a dump file for the mutation binary matrix.
     * @throws Exception
     */
    @Test
    public void generateSampleToGeneBinaryMatrix() throws Exception {
//        String dirName = "datasets/TCGA/GBM/FireHose/";
//        String mutationFileName = dirName + "FH_GBM.Mutation_Significance.Level_4.20120725/GBM.final_analysis_set.maf";
//        String clusterFileName = dirName + "2012FISubNetwork_Sample5_Modules_092012.txt";
        
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/OV.final_analysis_set.maf";
        String clusterFileName = dirName + "2009FISubNetwork_Sample3_Modules_092012.txt";
        
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationFileName, false);
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
        Set<String> cluster = clusters.get(3);
        List<String> clusterGenes = new ArrayList<String>(cluster);
        Collections.sort(clusterGenes);
        System.out.println("Total Genes: " + clusterGenes.size());
        System.out.println("Total Samples: " + sampleToGenes.size());
        
        // Split samples into two groups: one having mutated, and another no
        List<String> mutatedSamples = new ArrayList<String>();
        List<String> noMutatedSamples = new ArrayList<String>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(genes, clusterGenes);
            System.out.println(sample + "\t" + shared.size());
            if (shared.size() == 0)
                noMutatedSamples.add(sample);
            else
                mutatedSamples.add(sample);
        }
//        if (true)
//            return;
        Collections.sort(mutatedSamples);
        Collections.sort(noMutatedSamples);
        
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (String gene : clusterGenes)
            builder.append("\t").append(gene);
        builder.append("\n");
        for (String sample : mutatedSamples) {
            builder.append(sample);
            Set<String> genes = sampleToGenes.get(sample);
            for (String gene : clusterGenes) {
                if (genes.contains(gene))
                    builder.append("\t1");
                else
                    builder.append("\t0");
            }
            builder.append("\n");
        }
        for (String sample : noMutatedSamples) {
            builder.append(sample);
            Set<String> genes = sampleToGenes.get(sample);
            for (String gene : clusterGenes) {
                if (genes.contains(gene))
                    builder.append("\t1");
                else
                    builder.append("\t0");
            }
            builder.append("\n");
        }
        System.out.println(builder.toString());
    }
    
    public void loadTCGAFHDataSets(List<String> cancers,
                                    List<Map<String, Set<String>>> sampleToGenesList,
                                    List<String> clinFileNames) throws IOException {
        MATFileLoader matFileLoader = new MATFileLoader();
        //      TCGA OV
        cancers.add("OV");
        //     Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        //     String clinFileName = OVARIAN_DIR_NAME + "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt";
        String dirName = "datasets/TCGA/OvarianCancer/FireHose/";
        String mutationFileName = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/OV.final_analysis_set.maf";
        Map<String, Set<String>> sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        String clinFileName = dirName + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2012082500.0.0/OV.clin.transformed.txt";
        sampleToGenesList.add(sampleToGenes);
        clinFileNames.add(clinFileName);
        
        //      TCGA BRCA
        cancers.add("BRCA");
        dirName = "datasets/TCGA/BRCA/";
        mutationFileName = dirName + "FH_BRCA.Mutation_Significance.Level_4.20120725/BRCA.final_analysis_set.maf"; 
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        clinFileName = dirName + "FH_BRCA.Clinical_Pick_Tier1.Level_4.20120825/BRCA.clin.tranformed.txt";
        sampleToGenesList.add(sampleToGenes);
        clinFileNames.add(clinFileName);
        
        //      TCGA colorectal
        cancers.add("Colorectal");
        dirName = "datasets/TCGA/colorectal/";
        mutationFileName = dirName + "NonSynonymousSomaticMutations.txt";
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        clinFileName = dirName + "SimpleSurivalInfo.txt";
        sampleToGenesList.add(sampleToGenes);
        clinFileNames.add(clinFileName);
        
        // TCGA GBM
        cancers.add("GBM");
        dirName = "datasets/TCGA/GBM/FireHose/";
        mutationFileName = dirName + "FH_GBM.Mutation_Significance.Level_4.20120725/GBM.final_analysis_set.maf";
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        clinFileName = dirName + "FH_GBM.Clinical_Pick_Tier1.Level_4.20120825/GBM.clin.transformed.txt";
        sampleToGenesList.add(sampleToGenes);
        clinFileNames.add(clinFileName);
        
        //     // TCGA KIRC (Kidney)
        cancers.add("KIRC");
        dirName = "datasets/TCGA/KIRC/";
        mutationFileName = dirName + "FH_KIRC.Mutation_Significance.Level_4.20120725/KIRC.final_analysis_set.maf";
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        clinFileName = dirName + "FH_KIRC.Clinical_Pick_Tier1.Level_4.20120825/KIRC.clin.transformed.txt";
        sampleToGenesList.add(sampleToGenes);
        clinFileNames.add(clinFileName);
        
        // LAML (Leukemia)
        cancers.add("LAML");
        dirName = "datasets/TCGA/LAML/";
        mutationFileName = dirName + "gdac.broadinstitute.org_LAML.Mutation_Significance.Level_4.2012072500.0.0/LAML.final_analysis_set.maf";
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        sampleToGenesList.add(sampleToGenes);
        clinFileName = dirName + "gdac.broadinstitute.org_LAML.Clinical_Pick_Tier1.Level_4.2012082500.0.0/LAML.clin.transformed.120512.txt";
        clinFileNames.add(clinFileName);
        
        // UCEC (Uterine Corpus Endometrioid Carcinoma)
        cancers.add("UCEC");
        dirName = "datasets/TCGA/UCEC/";
        mutationFileName = dirName + "gdac.broadinstitute.org_UCEC.Mutation_Significance.Level_4.2012072500.0.0/UCEC.final_analysis_set.maf";
        sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        sampleToGenesList.add(sampleToGenes);
        clinFileName = dirName + "gdac.broadinstitute.org_UCEC.Clinical_Pick_Tier1.Level_4.2012082500.0.0/UCEC.clin.transformed.120512.txt";
        clinFileNames.add(clinFileName);
    }
    
    /**
     * This method is used to do MCL clustering using edge weights based on gene scores. Gene
     * scores can be generated from sample expression or MutSigs.
     * @throws Exception
     */
    @Test
    public void checkMCLClustering() throws Exception {
        // For survival test
        List<String> cancers = new ArrayList<String>();
        List<Map<String, Set<String>>> sampleToGeneList = new ArrayList<Map<String,Set<String>>>();
        List<String> clinFileNames = new ArrayList<String>();
        loadTCGAFHDataSets(cancers, sampleToGeneList, clinFileNames);
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        
        // MCL based on MutSig pvalues
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/";
//        String geneSigFile = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/" + 
//                "OV.sig_genes.txt";
//        Map<String, Double> geneToPvalue = loadGeneToMutSigPValue(geneSigFile);
//        Map<String, Double> geneToScore = new HashMap<String, Double>();
//        for (String gene : geneToPvalue.keySet()) {
//            Double pvalue = geneToPvalue.get(gene);
//            geneToScore.put(gene, -Math.log10(pvalue));
//        }
//        
        // Gene scores based on mutated sample numbers
        String mutationSampleFile = R3Constants.OVARIAN_DIR_NAME + "FireHose/" + 
                "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/" +
                "OV.final_analysis_set.maf";
        Map<String, Set<String>> sampleToGenes = new MATFileLoader().loadSampleToGenes(mutationSampleFile, 
                                                                                       false);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            geneToScore.put(gene, new Double(samples.size()));
        }
        
        MCLClusterWrapper mclHelper = new MCLClusterWrapper();
//        mclHelper.setInflation(7.5d);
        List<Set<String>> modules = mclHelper.mclClusterForGeneScores(geneToScore,
                                                                      false,
                                                                      null,
                                                                      null);
        for (int i = 0; i < modules.size(); i++) {
            Set<String> module = modules.get(i);
            if (module.size() < 3)
                continue; // Don't want to print out interactions.
            // Get average score
            Double total = 0.0d;
            for (String gene : module) {
                Double score = geneToScore.get(gene);
                if (score == null)
                    score = 0.0d;
                total += score;
            }
            total = total / module.size();
            System.out.println(i + "\t" + module.size() + "\t" + total);
            if (total > 6.0d) {
                System.out.println(module);
                doSurvivalAnalysisForMutationModule(survivalHelper, 
                                                    sampleToGeneList.get(0),
                                                    module, 
                                                    clinFileNames.get(0), 
                                                    "Module_" + i);
            }
        }
    }
    
}
