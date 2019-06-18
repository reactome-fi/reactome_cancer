/*
 * Created on Nov 6, 2009
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.commons.math.stat.inference.TestUtils;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.cluster.HierarchicalCluster;
import org.reactome.r3.cluster.HierarchicalClusterNode;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.graph.GraphComponent;
import org.reactome.r3.graph.GraphComponentSearchEngine;
import org.reactome.r3.graph.JGraphTUtilities;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.graph.TopologicalOverlapCalculator;
import org.reactome.r3.util.FisherExact;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

public class TCGAOvarianCancerClinicalAnalyzer extends TCGAOvarianCancerAnalyzer {
//    #private final String CLINICAL_FILE = OVARIAN_DIR_NAME + "data_110609/tcga_OV_clinical_csv_corrected.txt";
    private final String CLINICAL_FILE = OVARIAN_DIR_NAME + "data_031910/Batches9-22_tcga_OV_clinical_csv.2010-01-25-noDates.txt";
    //private final String CLINICAL_FILE = OVARIAN_DIR_NAME + "data_110609/ClinBatch17_19.txt";
    private final String TOTHILL_DIR = OVARIAN_DIR_NAME + "TothillDataset/";
    
    public TCGAOvarianCancerClinicalAnalyzer() {
    }
    
    /**
     * This method is used to check progress-free related genes from coxph filter in the 
     * genefilter package from R.
     * @throws Exception
     */
    @Test
    public void checkCoxphGenesFromR() throws Exception {
        String fileName = OVARIAN_DIR_NAME + "CoxphGenesFromGeneFilter.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        line = fu.readLine();
        // For TCGA
        String[] tokens = line.split(", ");
        List<String> tcgaGenes = Arrays.asList(tokens);
        line = fu.readLine();
        line = fu.readLine();
        List<String> tothillGenes = Arrays.asList(line.split(", "));
        System.out.println("TCGA genes: " + tcgaGenes.size());
        System.out.println("Tothill genes: " + tothillGenes.size());
        Set<String> shared = InteractionUtilities.getShared(tcgaGenes, tothillGenes);
        System.out.println("Shared genes: " + shared.size() + ": " + shared);
        // Check how many FI partners have been shared between these two gene sets
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fisInOVGenes = InteractionUtilities.getFIs(tcgaGenes, fis);
        System.out.println(fisInOVGenes);
        fu.saveInteractions(fisInOVGenes, 
                            OVARIAN_DIR_NAME + "FIsCoxphGenesFromGeneFilter.txt");
        if (true)
            return;
        shared = getSharedFIPartners(tcgaGenes, tothillGenes, fis);
        System.out.println("Shared partners: " + shared.size());
        // Want to do a random permutation test
        Map<String, Map<String, Double>> tcgaGeneExp = loadGeneExp();
        Set<String> tcgaAllGenes = new HashSet<String>(tcgaGeneExp.keySet());
        Map<String, Map<String, Double>> tothillGeneExp = loadTothillGeneExpData();
        Set<String> tothillAllGenes = new HashSet<String>(tothillGeneExp.keySet());
        int permutation = 10;
        List<Integer> permutationResults = new ArrayList<Integer>();
        Set<String> tcgaRandomGenes = null;
        Set<String> tothillRandomGenes = null;
        RandomData randomizer = new RandomDataImpl();
        for (int i = 0; i < permutation; i++) {
            tcgaRandomGenes = MathUtilities.randomSampling(tcgaAllGenes, 
                                                           tcgaGenes.size(),
                                                           randomizer);
            tothillRandomGenes = MathUtilities.randomSampling(tothillAllGenes,
                                                              tothillGenes.size(), 
                                                              randomizer);
            shared = getSharedFIPartners(tcgaRandomGenes,
                                         tothillRandomGenes,
                                         fis);
            System.out.println("Permutation " + i + ": " + shared.size());
            permutationResults.add(shared.size());
        }
        Collections.sort(permutationResults);
        for (Integer v : permutationResults)
            System.out.println(v);
    }
    
    private Set<String> getSharedFIPartners(Collection<String> genes1,
                                            Collection<String> genes2,
                                            Set<String> fis) {
        Set<String> tcgaFIs = InteractionUtilities.grepFIsContains(genes1, fis);
        Set<String> tcgaPartners = InteractionUtilities.grepIDsFromInteractions(tcgaFIs);
        tcgaPartners.removeAll(genes1);
        Set<String> tothillFIs = InteractionUtilities.grepFIsContains(genes2, fis);
        Set<String> tothillPartners = InteractionUtilities.grepIDsFromInteractions(tothillFIs);
        tothillPartners.removeAll(genes2);
        Set<String> shared = InteractionUtilities.getShared(tcgaPartners, tothillPartners);
        return shared;
    }
    
    /**
     * This method is used to generate a TOM for OV altered genes so that it can be loaded into R.
     * @throws Exception
     */
    @Test
    public void generateTOMForAlteredGenes() throws Exception {
        Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        System.out.println("Total genes: " + geneToSamples.size());
        List<String> geneList = new ArrayList<String>(geneToSamples.keySet());
        Collections.sort(geneList);
        
        // Load FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        geneList.retainAll(fiGenes);
        System.out.println("Total genes in FIs: " + geneList.size());
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        
        // Prepare topological overlap calculation
        TopologicalOverlapCalculator topCalculator = new TopologicalOverlapCalculator();
        topCalculator.setIdToPartners(geneToPartners);
        
        // Output the matrix
        String outputFileName = OVARIAN_DIR_NAME + "TOMForAlteredGenes011110.txt";
        fu.setOutput(outputFileName);
        // Generate gene list
        StringBuilder line = new StringBuilder();
        for (String gene : geneList)
            line.append("\t").append(gene);
        fu.printLine(line.toString());
        line.setLength(0);
        for (String gene : geneList) {
            line.append(gene);
            for (String gene1 : geneList) {
                double dist = topCalculator.calculateDistance(gene, gene1);
                line.append("\t").append(dist);
            }
            fu.printLine(line.toString());
            line.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This method is used to hierarchically cluster a set of genes.
     * @throws Exception
     */
    @Test
    public void hierarchicalClusterCorrelatedGenes() throws Exception {
        List<String> correlatedGenes = loadPlatCorrelatedGenes();
        System.out.println("Total correlated genes: " + correlatedGenes.size());
        hierarchicalClusterGenes(correlatedGenes);
    }
    
    @Test
    public void compareAlteredGenesInPlatDiffSamples() throws Exception {
        Map<String, Double> sampleToPlatFreeInt = loadSampleToPlatinumFreeInterval(true);
        System.out.println("Total samples: " + sampleToPlatFreeInt.size());
        Set<String> platResSamples = new HashSet<String>();
        Set<String> platSenSamples = new HashSet<String>();
        splitPlatDiffSamples(sampleToPlatFreeInt, 
                             platResSamples,
                             platSenSamples);
        
        // Used to get shared samples
        //Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(true);
        Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        
        Set<String> genesInPlatResSamples = new HashSet<String>();
        Set<String> genesInPlatSenSamples = new HashSet<String>();
        for (String sample : platResSamples) {
            Set<String> genes = sampleToGenes.get(sample);
            if (genes != null)
                genesInPlatResSamples.addAll(genes);
        }
        for (String sample : platSenSamples) {
            Set<String> genes = sampleToGenes.get(sample);
            if (genes != null)
                genesInPlatSenSamples.addAll(genes);
        }
        System.out.println("Total genes in plat res samples: " + genesInPlatResSamples.size());
        for (String gene : genesInPlatResSamples)
            System.out.println(gene);
        System.out.println("Total genes in plat sen samples: " + genesInPlatSenSamples.size());
        for (String gene : genesInPlatSenSamples)
            System.out.println(gene);
        Set<String> shared = InteractionUtilities.getShared(genesInPlatResSamples, genesInPlatSenSamples);
        System.out.println("Shared genes: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(geneToSamples.size(),
                                                                    genesInPlatSenSamples.size(), 
                                                                    genesInPlatResSamples.size(), 
                                                                    shared.size());
        System.out.println("Pvalue for sharing: " + pvalue);
        FisherExact fisher = new FisherExact(geneToSamples.size());
        pvalue = fisher.getLeftTailedP(shared.size(),
                                       genesInPlatResSamples.size() - shared.size(),
                                       genesInPlatResSamples.size() - shared.size(),
                                       geneToSamples.size() - genesInPlatResSamples.size() - genesInPlatResSamples.size() + shared.size());
        System.out.println("Left tailed fisher exact: " + pvalue);
        pvalue = fisher.getRightTailedP(shared.size(),
                                       genesInPlatResSamples.size() - shared.size(),
                                       genesInPlatResSamples.size() - shared.size(),
                                       geneToSamples.size() - genesInPlatResSamples.size() - genesInPlatResSamples.size() + shared.size());
        System.out.println("Right tailed fisher exact: " + pvalue);
    }

    private void splitPlatDiffSamples(Map<String, Double> sampleToPlatFreeInt,
                                      Set<String> platResSamples,
                                      Set<String> platSenSamples) {
        for (String sample : sampleToPlatFreeInt.keySet()) {
            Double value = sampleToPlatFreeInt.get(sample);
            if (value <= 6)
                platResSamples.add(sample);
            else if (value >= 9)
                platSenSamples.add(sample);
        }
        System.out.println("Total plat sen: " + platSenSamples.size());
        System.out.println("Total plat res: " + platResSamples.size());
    }
    
    private void divideSamplesIntoUpperLowerOnSurvival(Set<String> longSamples,
                                                       Set<String> shortSamples,
                                                       Map<String, Set<String>> sampleToGenes) throws IOException {
        Map<String, Double> sampleToRate = loadSampleToSurvivalRate(true);
        // Get shared samples
        Set<String> sharedSamples = new HashSet<String>(sampleToGenes.keySet());
        sharedSamples.retainAll(sampleToRate.keySet());
        double lower = 33.80d;
        double upper = 33.80d;
        for (String sample : sharedSamples) {
            Double value = sampleToRate.get(sample);
            if (value < lower)
                shortSamples.add(sample);
            else if (value > upper)
                longSamples.add(sample);
        }
        System.out.println("Long samples: " + longSamples.size());
        System.out.println("Short samples: " + shortSamples.size());
        System.out.println("Total samples: " + sharedSamples.size());
    }
    
    /**
     * This method is used to check a single module and its genes for distinguishing sample survival rates.
     * @throws Exception
     */
    @Test
    public void calculateRelationForAlteredGenesOnSurvivalRate() throws Exception {
        Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        filterSampleToGenes(sampleToGenes, 2);
        Set<String> longSamples = new HashSet<String>();
        Set<String> shortSamples = new HashSet<String>();
        divideSamplesIntoUpperLowerOnSurvival(longSamples, 
                                              shortSamples, 
                                              sampleToGenes);
        
        // Check based on modules
        //String clusterFileName = OVARIAN_DIR_NAME  + "ClustersForMutatedGenesSpectral_2Sample_040110.txt";
        //String clusterFileName = OVARIAN_DIR_NAME + "ClustersForMutatedGenesEdge_2Sample_040210.txt";
        //NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        //List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
//        filterSampleToGenes(sampleToGenes, 2);
        Set<String> allGenes1 = InteractionUtilities.grepAllGenes(sampleToGenes);
        List<Set<String>> clusters = new SpectralPartitionNetworkCluster().clusterGenes(allGenes1);
        
        System.out.println("Module\tSamplesInUpper\tSamplesInLower\tTwo_tail_P\tLeft_P\tRight_P");
        int a, b;
        FisherExact fisherTest = new FisherExact(sampleToGenes.size());
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            // Check samples in cluster 2
            
            a = b = 0;
            for (String sample : longSamples) {
                Set<String> genes = sampleToGenes.get(sample);
                Set<String> shared = InteractionUtilities.getShared(genes, cluster);
                if (shared.size() > 0)
                    a ++;
            }
            for (String sample : shortSamples) {
                Set<String> genes = sampleToGenes.get(sample);
                Set<String> shared = InteractionUtilities.getShared(genes, cluster);
                if (shared.size() > 0)
                    b ++;
            }
            double pvalue = fisherTest.getTwoTailedP(a, b, longSamples.size() - a, shortSamples.size() - b); 
            double left = fisherTest.getLeftTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            double right = fisherTest.getRightTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            
            System.out.println("Module " + i + "\t" + a + "\t" + b + "\t" + String.format("%.3f", pvalue) + "\t" + 
                               String.format("%.3f", left) + "\t" + String.format("%.3f", right));
        }
        // Check with a network modulation
        List<Integer> checkingModules = new ArrayList<Integer>();
        //checkingModules.add(1);
        checkingModules.add(2);
        checkingModules.add(9);
        a = b = 0;
        Set<String> checkingGenes = new HashSet<String>();
        for (Integer index : checkingModules)
            checkingGenes.addAll(clusters.get(index));
        for (String sample : longSamples) {
            Set<String> genes = sampleToGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(genes, checkingGenes);
            if (shared.size() > 0)
                a ++;
        }
        for (String sample : shortSamples) {
            Set<String> genes = sampleToGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(genes, checkingGenes);
            if (shared.size() > 0)
                b ++;
        }
        double pvalue = fisherTest.getTwoTailedP(a, b, longSamples.size() - a, shortSamples.size() - b); 
        double left = fisherTest.getLeftTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
        double right = fisherTest.getRightTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
        
        System.out.println("Module " + checkingModules + "\t" + a + "\t" + b + "\t" + String.format("%.3f", pvalue) + "\t" + 
                           String.format("%.3f", left) + "\t" + String.format("%.3f", right));
        if (true)
            return;
        
//        List<String> geneList = new ArrayList<String>(cluster);
        System.out.println("\nGene\tSamplesInUpper\tSamplesInLower\tTwo_tail_P\tLeft_P\tRight_P");
        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
        List<String> geneList = new ArrayList<String>(allGenes);
        Collections.sort(geneList);
        for (String gene : geneList) {
            a = b = 0;
            for (String sample : longSamples) {
                Set<String> genes = sampleToGenes.get(sample);
                if (genes.contains(gene))
                    a ++;
            }
            for (String sample : shortSamples) {
                Set<String> genes = sampleToGenes.get(sample);
                if (genes.contains(gene))
                    b ++;
            }
            pvalue = fisherTest.getTwoTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            left = fisherTest.getLeftTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            right = fisherTest.getRightTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            System.out.println(gene + "\t" + a + "\t" + b + "\t" + String.format("%.3f", pvalue) + "\t" + 
                               String.format("%.3f", left) + "\t" + String.format("%.3f", right));
        }
//        // Do a t-test
//        List<Double> mutatedRates = new ArrayList<Double>();
//        List<Double> nonMutatedRates = new ArrayList<Double>();
//        for (String sample : sharedSamples) {
//            Set<String> genes = sampleToGenes.get(sample);
//            Double rate = sampleToRate.get(sample);
//            Set<String> shared = InteractionUtilities.getShared(genes, cluster);
//            if (shared.size() > 0)
//                mutatedRates.add(rate);
//            else
//                nonMutatedRates.add(rate);
//        }
//        System.out.println("Total sample having mutated genes in cluster: " + mutatedRates.size());
//        System.out.println(mutatedRates);
//        System.out.println("Total sample having no mutated genes in cluster: " + nonMutatedRates.size());
//        System.out.println(nonMutatedRates);
//        pvalue = MathUtilities.calculateTTest(mutatedRates, nonMutatedRates);
//        System.out.println("pvalue from ttest: " + pvalue);
    }
    
    /**
     * This method is used to check the relation between individual sequence altered genes
     * and plat types. The calculation is based on Fisher exact test.
     * @throws Exception
     */
    @Test
    public void calculateRelationForAlteredGenesOnPlatInfo() throws Exception {
        Map<String, Double> sampleToPlatFreeInt = loadSampleToPlatinumFreeInterval(true);
        System.out.println("Total samples: " + sampleToPlatFreeInt.size());
        
        // Used to get shared samples
        //Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(true);
        Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        // The following statements are used to ensure only shared samples are used in two maps.
        Set<String> keys = sampleToGenes.keySet();
        keys.retainAll(sampleToPlatFreeInt.keySet());
        sampleToPlatFreeInt.keySet().retainAll(keys);
        
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        
        Set<String> platResSamples = new HashSet<String>();
        Set<String> platSenSamples = new HashSet<String>();
        splitPlatDiffSamples(sampleToPlatFreeInt, 
                             platResSamples,
                             platSenSamples);
        Set<String> allSamples = new HashSet<String>();
        allSamples.addAll(platSenSamples);
        allSamples.addAll(platResSamples);
        sampleToGenes.keySet().retainAll(allSamples);
        
        int total = platSenSamples.size() + platResSamples.size();
        //System.out.println("\nGene\tPvalueForPlatSen\tPvalueForPlatRes");
        Set<String> selectedGenes = new HashSet<String>();
        double pvalueCutoff = 0.25;
        String fileName = OVARIAN_DIR_NAME + "PlatSignificantAlteredGenesPValue_homo_mutated_050.arff";
        Set<String> platSenGenes = new HashSet<String>();
        Set<String> platResGenes = new HashSet<String>();
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            Set<String> platSenSamplesInGene = InteractionUtilities.getShared(samples, platSenSamples);
            Set<String> platResSamplesInGene = InteractionUtilities.getShared(samples, platResSamples);
            // For platinum resistant
            double pvalue1 = MathUtilities.calculateHypergeometricPValue(total, 
                                                                         samples.size(),
                                                                         platSenSamples.size(), 
                                                                         platSenSamplesInGene.size());
            if (pvalue1 <= pvalueCutoff) {
                selectedGenes.add(gene);
                platSenGenes.add(gene);
            }
            // For platinum sensitive
                double pvalue2 = MathUtilities.calculateHypergeometricPValue(total, 
                                                                             samples.size(), 
                                                                             platResSamples.size(),
                                                                             platResSamplesInGene.size());
            if (pvalue2 <= pvalueCutoff) {
                selectedGenes.add(gene);
                platResGenes.add(gene);
            }
            //System.out.println(gene + "\t" + pvalue1 + "\t" + pvalue2);
        }
        
//        // Compare genes from graph component
//        String fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatPos121809.txt";
//        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
//        List<GraphComponent> components = componentAnalyzer.loadComponentsInOrder(fileName);
//        System.out.println("\nCheck genes shared with components:");
//        PathwayAnnotator annotator = new PathwayAnnotator();
//        annotator.setFDRThreshold(0.50);
//        //annotator.setPValueThreshold(0.001);
//        for (int i = 0; i < components.size(); i++) {
//            GraphComponent comp = components.get(i);
//            Set<String> shared1 = InteractionUtilities.getShared(comp.getAllNodes(),
//                                                                 selectedGenes);
//            double pvalue = MathUtilities.calculateHypergeometricPValue(10000, 
//                                                                        comp.getAllNodes().size(),
//                                                                        selectedGenes.size(), 
//                                                                        shared1.size());
//            System.out.println("Shared genes between component " + i + ": " + shared1.size() + ", " + pvalue);
//            //annotator.annotateGenesWithFDR(comp.getAllNodes());
//            if (i == 2)
//                break;
//        }
        
//        // Generate an arff file
//        fu.setInput(OVARIAN_DIR_NAME + "PlatSignificantAlteredGenes.txt");
//        List<String> selectedGenes = new ArrayList<String>();
//        String line = null;
//        while ((line = fu.readLine()) != null)
//            selectedGenes.add(line);
//        fu.close();

        System.out.println("\nTotal selected genes: " + selectedGenes.size());
//        for (String gene : selectedGenes)
//            System.out.println(gene);
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.20);
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            genes.retainAll(selectedGenes);
        }
//        annotator.annotatePathwayBasedOnSamples(sampleToGenes, 
//                                                false,
//                                                true);
//        System.out.println("\nTotal plat sen genes: " + platSenGenes.size());
//        Map<String, Set<String>> copy = new HashMap<String, Set<String>>(sampleToGenes);
//        copy.keySet().retainAll(platSenSamples);
//        System.out.println("Total plat sen samples for checking: " + copy.size());
//        annotator.annotatePathwayBasedOnSamples(copy,
//                                                false,
//                                                true);
//        System.out.println("\nTotal plat res genes: " + platResGenes.size());
        Map<String, Set<String>> copy = new HashMap<String, Set<String>>(sampleToGenes);
        copy.keySet().retainAll(platResSamples);
        System.out.println("Total plat res samples for checking: " + copy.size());
//        annotator.annotatePathwayBasedOnSamples(copy,
//                                                false,
//                                                true);
        if (true)
            return;

        List<String> geneList = new ArrayList<String>(selectedGenes);
        Collections.sort(geneList);
        generateARFFForAlteredGenes(sampleToGenes, 
                                    platResSamples,
                                    platSenSamples,
                                    geneList, 
                                    fileName);
    }

    @Test
    public void generateTestSampleToAlteredGenes() throws Exception {
        TCGAOvarianCancerClinicalAnalyzer clinicalAnalyzer = new TCGAOvarianCancerClinicalAnalyzer();
        Map<String, Double> sampleToPlatFreeInt = clinicalAnalyzer.loadSampleToPlatinumFreeInterval(true);
        Map<String, Double> fullSampleToPlatFreeInt = clinicalAnalyzer.loadSampleToPlatinumFreeInterval(false);
        fullSampleToPlatFreeInt.keySet().removeAll(sampleToPlatFreeInt.keySet());
        sampleToPlatFreeInt = fullSampleToPlatFreeInt;
        System.out.println("Total samples: " + sampleToPlatFreeInt.size());
        
        // Used to get shared samples
        //Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        // The following statements are used to ensure only shared samples are used in two maps.
        Set<String> keys = sampleToGenes.keySet();
        keys.retainAll(sampleToPlatFreeInt.keySet());
        sampleToPlatFreeInt.keySet().retainAll(keys);
        
        Set<String> platResSamples = new HashSet<String>();
        Set<String> platSenSamples = new HashSet<String>();
        splitPlatDiffSamples(sampleToPlatFreeInt, platResSamples,
                             platSenSamples);
        Set<String> selectedGenes = fu.loadInteractions(OVARIAN_DIR_NAME + "PlatSignificantAlteredGenesPValue_025.txt");
        String fileName = OVARIAN_DIR_NAME + "PlatSignificantAlteredGenesPValue_025_followup.arff";
        List<String> geneList = new ArrayList<String>(selectedGenes);
        Collections.sort(geneList);
        generateARFFForAlteredGenes(sampleToGenes, 
                                    platResSamples, 
                                    platSenSamples, 
                                    geneList, 
                                    fileName);
    }
    
    @Test
    public void generateFIsForGSEA() throws Exception {
        Map<String, Double> sampleToPlat = loadSampleToPlatinumFreeInterval(true);
        // Divide the samples into two parts
        List<String> platResSamples = new ArrayList<String>();
        List<String> platSenSamples = new ArrayList<String>();
        for (String sample : sampleToPlat.keySet()) {
            Double plat = sampleToPlat.get(sample);
            if (plat <= 6.0)
                platResSamples.add(sample);
            else if (plat >= 9.0)
                platSenSamples.add(sample);
        }
        // Generate sample file
        String outFileName = OVARIAN_DIR_NAME + "PlatSamples_Batch9_15.cls";
        fu.setOutput(outFileName);
        fu.printLine((platResSamples.size() + platSenSamples.size()) + " 2 1");
        fu.printLine("# RES SEN");
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < platResSamples.size(); i++)
            builder.append(0).append(" ");
        for (int i = 0; i < platSenSamples.size(); i++)
            builder.append(1).append(" ");
        fu.printLine(builder.toString());
        fu.close();
        Map<String, Map<String, Double>> geneExp = loadGeneExp();
        System.out.println("Total genes in exp: " + geneExp.size());
        // Try to find the smallest p-value
        outFileName = OVARIAN_DIR_NAME + "PlatSamplesGSEAExp_Batch9_15.txt";
        // Generate the title
        fu.setOutput(outFileName);
        builder.setLength(0);
        builder.append("Name\tDESCRIPTION");
        for (String sample : platResSamples)
            builder.append("\t").append(sample);
        for (String sample : platSenSamples)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        for (String gene : geneExp.keySet()) {
            builder.setLength(0);
            builder.append(gene).append("\tNA");
            Map<String, Double> sampleToExp = geneExp.get(gene);
            for (int i = 0; i < platResSamples.size(); i++) {
                Double v = sampleToExp.get(platResSamples.get(i));
                builder.append("\t").append(v);
            }
            for (int i = 0; i < platSenSamples.size(); i++) {
                Double v = sampleToExp.get(platSenSamples.get(i));
                builder.append("\t").append(v);
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * This method is used to check FIs for gene expression data set. The method runs very slow!
     * @throws Exception
     */
    @Test
    public void checkFIsInGeneExp() throws Exception {
        Map<String, Double> sampleToPlat = loadSampleToPlatinumFreeInterval(true);
        // Divide the samples into two parts
        List<String> platResSamples = new ArrayList<String>();
        List<String> platSenSamples = new ArrayList<String>();
        for (String sample : sampleToPlat.keySet()) {
            Double plat = sampleToPlat.get(sample);
            if (plat <= 6.0)
                platResSamples.add(sample);
            else if (plat >= 9.0)
                platSenSamples.add(sample);
        }
        double[] resValues = new double[platResSamples.size()];
        double[] senValues = new double[platSenSamples.size()];
        //Map<String, Map<String, Double>> geneExp = loadGeneExp();
        String expFileName = "datasets/TCGA/OvarianCancer/GeneExp090209/broad.mit.edu__HT_HG-U133A__gene_expression_analysis";
        TCGADataHelper dataHelper = new TCGADataHelper();
        Map<String, Map<String, Double>> geneExp = dataHelper.loadGeneExpFromPortal(expFileName);
        dataHelper.zTransformBasedOnSample(geneExp);
        System.out.println("Total genes in exp: " + geneExp.size());
        // Try to find the smallest p-value
        double minP = Double.MAX_VALUE;
        double pvalue;
        long time1 = System.currentTimeMillis();
        for (String gene : geneExp.keySet()) {
            Map<String, Double> sampleToExp = geneExp.get(gene);
            for (int i = 0; i < platResSamples.size(); i++) {
                Double v = sampleToExp.get(platResSamples.get(i));
                resValues[i] = v;
            }
            for (int i = 0; i < platSenSamples.size(); i++) {
                Double v = sampleToExp.get(platSenSamples.get(i));
                senValues[i] = v;
            }
            pvalue = TestUtils.tTest(resValues, senValues);
            if (pvalue < minP)
                minP = pvalue;
        }
        System.out.println("Min pvalue for gene: " + minP);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for genes ttest: " + (time2 - time1));
        if (true)
            return;
        //Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        //Set<String> fisInGenes = InteractionUtilities.getFIs(geneExp.keySet(), fis);
        // Save in a file to run it again
        String geneFIFileName = OVARIAN_DIR_NAME + "FIsInGeneExpBatch9_15.txt";
        //fu.saveInteractions(fisInGenes, geneFIFileName);
        Set<String> fisInGenes = fu.loadInteractions(geneFIFileName);
        System.out.println("Total FIs in exp: " + fisInGenes.size());
        int index = 0;
        String gene1, gene2, sample;
        Map<String, Double> sampleToExp1, sampleToExp2;
        minP = Double.MAX_VALUE;
        String outFileName = OVARIAN_DIR_NAME + "PValuesFromTForFIsInGeneExpBatch9_15.txt";
        fu.setOutput(outFileName);
        fu.printLine("FI\tP-value");
        for (String fi : fisInGenes) {
            index = fi.indexOf("\t");
            gene1 = fi.substring(0, index);
            sampleToExp1 = geneExp.get(gene1);
            gene2 = fi.substring(index + 1);
            sampleToExp2 = geneExp.get(gene2);
            // Get average for both values
            for (int i = 0; i < platResSamples.size(); i++) {
                sample = platResSamples.get(i);
                Double v1 = sampleToExp1.get(sample);
                Double v2 = sampleToExp2.get(sample);
                resValues[i] = (v1 + v2) / 2.0;
            }
            for (int i = 0; i < platSenSamples.size(); i++) {
                sample = platSenSamples.get(i);
                Double v1 = sampleToExp1.get(sample);
                Double v2 = sampleToExp2.get(sample);
                senValues[i] = (v1 + v2) / 2.0;
            }
            pvalue = TestUtils.tTest(resValues, senValues);
            if (pvalue < minP)
                minP = pvalue;
            fu.printLine(gene1 + ":" + gene2 + "\t" + pvalue);
        }
        fu.close();
        System.out.println("Min pvalue for FIs: " + minP);
        long time3 = System.currentTimeMillis();
        System.out.println("Total time for FIs: " + (time3 - time2));
    }
    
    private void generateARFFForAlteredGenes(Map<String, Set<String>> sampleToGenes,
                                             Set<String> platResSamples,
                                             Set<String> platSenSamples,
                                             Collection<String> selectedGenes,
                                             String fileName) throws IOException {
        fu.setOutput(fileName);
        fu.printLine("@relation SingleGeneBasedPlat");
        fu.printLine("");
        fu.printLine("@attribute PlatinumResistance {true, false}");
        for (String gene : selectedGenes)
            fu.printLine("@attribute " + gene + " {true, false}");
        fu.printLine("");
        fu.printLine("@data");
        StringBuilder builder = new StringBuilder();
        for (String sample : platResSamples) {
            builder.append(true);
            Set<String> genes = sampleToGenes.get(sample);
            for (String gene : selectedGenes) {
                builder.append(",").append(genes.contains(gene));
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        for (String sample : platSenSamples) {
            builder.append(false);
            Set<String> genes = sampleToGenes.get(sample);
            for (String gene : selectedGenes) {
                builder.append(",").append(genes.contains(gene));
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This method is used to check correlation for clusters from hierarchical clustering.
     */
    @Test
    public void checkCorrelationForHierarchicalClusters() throws Exception {
        HierarchicalCluster clustering = new HierarchicalCluster();
        List<HierarchicalClusterNode> clusters = clustering.loadHierarchicalClusters(OVARIAN_DIR_NAME + "HierarchicalClusterOfPlatCorGenes.txt");
        // Load these data sets
        Map<String, Double> sampleToPlatInterval = loadSampleToPlatinumFreeInterval();
        Map<String, Map<String, Double>> geneToExp = loadGeneExp();
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        scoreCalculator.setGeneToExp(geneToExp);
        scoreCalculator.setSampleToPlatinum(sampleToPlatInterval);
        for (HierarchicalClusterNode cluster : clusters) {
            Set<String> genes = cluster.getIds();
            if (genes.size() < 5)
                continue;
            Map<String, Double> values = scoreCalculator.calculateAverage(genes);
            Double score = scoreCalculator.calculateCorBetweenExpAndPlatInverval(values, 
                                                                                 sampleToPlatInterval);
            if (Math.abs(score) < 0.15)
                continue;
            System.out.println(cluster.pathDistance + ", " + 
                               cluster.ids.size() + ", " + cluster.ids + "\n" +
                               score);
            // Want to print out score for each genes
            for (String gene : genes) {
                values = geneToExp.get(gene);
                score = scoreCalculator.calculateCorBetweenExpAndPlatInverval(values,
                                                                              sampleToPlatInterval);
                System.out.println("\t" + gene + ": " + score);
            }
        }
        String outFileName = OVARIAN_DIR_NAME + "SampleExpToPlatIntervalHierPosNegCluster.arff";
        generateARFFForHiClusters(scoreCalculator,
                                  outFileName);
    }

    private void generateARFFForHiClusters(PlatinumFreeIntervalScoreCalculator scoreCalculator,
                                           String fileName) throws IOException {
        // The following two clusters are extracted from hiearchical clustering with modification
        // The last three genes in the first negative cluster were picked from the second positive
        // cluster since they have negative correlation
        String negativeCluster = "EXOC2, RHOD, RHOT1, PIK3R3, PIP5K1A, ARF5, EPB41L1, ZYX, EFNA1";
        String positiveCluster = "SHC3, CD79A, LMO4, PRKCE, SPRY2, SPRY4, TDGF1, SMURF2, SPRY1";
        List<String> list = new ArrayList<String>();
        list.add(negativeCluster);
        list.add(positiveCluster);
        Map<String, Double> negSampleToGeneExp = null;
        Map<String, Double> posSampleToGeneExp = null;
        Map<String, Double> sampleToPlatInterval = scoreCalculator.getSampleToPlatinum();
        for (String cluster : list) {
            String[] tokens = cluster.split(", ");
            List<String> genes = Arrays.asList(tokens);
            Map<String, Double> values = scoreCalculator.calculateAverage(genes);
            if (cluster.equals(negativeCluster))
                negSampleToGeneExp = values;
            else
                posSampleToGeneExp = values;
            Double score = scoreCalculator.calculateCorBetweenExpAndPlatInverval(values, 
                                                                                 sampleToPlatInterval);
            System.out.println(cluster + ": " + score);
            //printGeneExpToPlatInterval(sampleToPlatInterval, 
            //                           values);
        }
        System.out.println("\nCombined output from positive and negative hierarchical clusters in the arff format:");
        List<String> attNames = new ArrayList<String>();
        attNames.add("PosGeneExp");
        attNames.add("NegGeneExp");
        List<Map<String, Double>> expData = new ArrayList<Map<String,Double>>();
        expData.add(posSampleToGeneExp);
        expData.add(negSampleToGeneExp);
        generateARFF(scoreCalculator,
                    expData,
                    attNames,
                    fileName);
    }
    
    /**
     * Generate an arff file for gene expression and platimum.
     * @param sampleToPlatInterval
     * @param geneExpData should have the same order as attNames.
     * @param attNames
     * @param fileName
     */
    private void generateARFF(PlatinumFreeIntervalScoreCalculator scoreCalcuator,
                              List<Map<String, Double>> geneExpData,
                              List<String> attNames,
                              String fileName) throws IOException {
        fu.setOutput(fileName);
        printArffHeader(attNames, false);
        StringBuilder builder = new StringBuilder();
        Map<String, Double> sampleToPlatInterval = scoreCalcuator.getSampleToPlatinum();
        for (String sample : sampleToPlatInterval.keySet()) {
            // A quick check to make sure this sample has expression value
            Map<String, Double> tmpExp = geneExpData.get(0);
            if (!scoreCalcuator.hasExpressionValue(sample, tmpExp)) {
                continue;
            }
            builder.setLength(0);
            Double platInt = sampleToPlatInterval.get(sample);
            if (platInt <= 6.0) {
                builder.append(true);
            }
            else if (platInt >= 9.0)
                builder.append(false);
            else
                continue; // Don't want to output samples in between.
            for (Map<String, Double> sampleToExp : geneExpData) {
                Double exp = scoreCalcuator.getExpressionValue(sample, sampleToExp);
                builder.append(",").append(exp);
            }
            fu.printLine(builder.toString());
            System.out.println(sample + "\t" + sampleToPlatInterval.get(sample));
        }
        fu.close();
    }

    private void printArffHeader(List<String> attNames,
                                 boolean isBoolean) throws IOException {
        fu.printLine("@relation Ovarian_Cancer_Platinum_Free_Interval\n");
        fu.printLine("@attribute PlatinumResistence {true, false}");
        //fu.printLine("@attribute PlatinumResistence numeric");
        for (String attName : attNames) {
            // Always should be numeric
            if (isBoolean)
                fu.printLine("@attribute " + attName + " {true, false}");
            else
                fu.printLine("@attribute " + attName + " numeric");
        }
        fu.printLine("\n@data");
    }
    
    /**
     * A method is used to search a graph component.
     * @throws Exception
     */
    @Test
    public void searchForGraphComponents() throws Exception {
        Map<String, Map<String, Double>> geneToExpData = loadGeneExp();
        Map<String, Double> sampleToPlatinum = loadSampleToPlatinumFreeInterval(true);
        // Want to work on genes with expression data only
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> geneInFis = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> keyset = geneToExpData.keySet();
        keyset.retainAll(geneInFis);
        Set<String> fisInExp = InteractionUtilities.getFIs(keyset, fis);
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        Graph<String, DefaultEdge> graph = JGraphTUtilities.createGraph(keyset, fisInExp);
        // For search
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        // Try to find negative component
        scoreCalculator.setIsForNegative(false);
        String outFileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatPos121809.txt";
        scoreCalculator.setGeneToExp(geneToExpData);
        scoreCalculator.setSampleToPlatinum(sampleToPlatinum);
        GraphComponentSearchEngine searcher = new GraphComponentSearchEngine();
        searcher.setScoreCalculator(scoreCalculator);
        searcher.setMaxDepth(2);
        searcher.setSearchDepth(1);
        searcher.search(graph);
        List<GraphComponent> components = searcher.getFoundComponents();
        System.out.println("Total componenents: " + components.size());
        //String outFileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlat.txt";
        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
        componentAnalyzer.outputFoundComponents(outFileName,
                                                components);
    }
    
    /**
     * Check some features in the first graph component.
     * @throws Exception
     */
    @Test
    public void checkCorrelationForGraphComponent() throws Exception {
        //PlatinumFreeIntervalScoreCalculator scoreCalculator = initScoreCalculator();
        
        String fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatPos121809.txt";
        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
        List<GraphComponent> components = componentAnalyzer.loadComponentsInOrder(fileName);
        // Check gene overlapes with the first two components
        GraphComponent comp1 = components.get(0);
        GraphComponent comp2 = components.get(1);
        Set<String> shared = InteractionUtilities.getShared(comp1.getAllNodes(), 
                                                            comp2.getAllNodes());
        System.out.println("Component 0: " + comp1.getAllNodes().size());
        System.out.println("Component 1: " + comp2.getAllNodes().size());
        System.out.println("Shared genes in the first two components: " + shared.size());
//        for (int i = 0; i < components.size(); i++) {
//            // Just want to check the first component
//            GraphComponent first = components.get(i);
//            double score = scoreCalculator.calculateScore(first);
//            System.out.println("Component " + i + "\t" + score + "\t" + first.getAllNodes().size());
//            if (i == 10)
//                break;
////            // Used to create scatter plot
////            Map<String, Double> firstSampleToExp = scoreCalculator.calculateAverage(first.getAllNodes());
////            printGeneExpToPlatInterval(sampleToPlatinum,
////                                       firstSampleToExp);
//        }
    }
    
    /**
     * This method is used to generate an arff file for GraphComponent found using a search algorithm.
     * @throws Exception
     */
    @Test
    public void generateARFFFileForGraphComponents() throws Exception {
        PlatinumFreeIntervalScoreCalculator scoreCalculator = initScoreCalculator();
        
        String fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlat.txt";
        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
        List<GraphComponent> components = componentAnalyzer.loadComponentsInOrder(fileName);
        List<Map<String, Double>> componentExpData = new ArrayList<Map<String,Double>>();
        List<String> attNames = new ArrayList<String>();
        // Try the first two graph components.
        for (int i = 0; i < 0; i++) {
            GraphComponent comp = components.get(i);
            Set<String> genes = comp.getAllNodes();
            Map<String, Double> compExp = scoreCalculator.calculateAverage(genes);
            componentExpData.add(compExp);
            attNames.add("PosComponent" + i);
            // Just a quick check
            double score = scoreCalculator.calculateScore(comp);
            System.out.println("Score for component" + i + ": " + score);
        }
        // For negative part
        fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatNeg.txt";
        components = componentAnalyzer.loadComponentsInOrder(fileName);
        // Try the first two graph components.
        for (int i = 0; i < 2; i++) {
            GraphComponent comp = components.get(i);
            Set<String> genes = comp.getAllNodes();
            Map<String, Double> compExp = scoreCalculator.calculateAverage(genes);
            componentExpData.add(compExp);
            attNames.add("NegComponent" + i);
            // Just a quick check
            double score = scoreCalculator.calculateScore(comp);
            System.out.println("Score for component" + i + ": " + score);
        }
        String outputFileName = OVARIAN_DIR_NAME + "SampleExpToPlatIntervalGraphComponentsNeg.arff";
        generateARFF(scoreCalculator,
                     componentExpData,
                     attNames,
                     outputFileName);
    }
    
    @Test
    public void generateARFFFileForNetworkModules() throws Exception {
        Map<String, Double> sampleToPlatInterval = loadSampleToPlatinumFreeInterval(true);
        Map<String, Double> fullSampleToPlatInterval = loadSampleToPlatinumFreeInterval(false);
        fullSampleToPlatInterval.keySet().removeAll(sampleToPlatInterval.keySet());
        sampleToPlatInterval = fullSampleToPlatInterval;
        // This file was generated from TCGAOvarianCancerAnalyzer.generateSampleToNetworkModuleMatrixInDegrees()
        //String fileName = OVARIAN_DIR_NAME + "SampleToTOMNetworkModulesForAlteredGenes011110.txt";
        //String fileName = OVARIAN_DIR_NAME + "SampleToMCLNetworkModulesForAlteredGenes112309.txt";
        String fileName = OVARIAN_DIR_NAME  + "SamplesToClustersForMutatedGenesSpectral_2Sample_040110.txt";
        boolean isBoolean = true;
        Map<String, List<Object>> sampleToValues = new HashMap<String, List<Object>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        List<String> headerList = new ArrayList<String>();
        for (int i = 1; i < headers.length; i++) {
            String header = headers[i];
            header = header.replaceAll("( )|\\(|\\)|\\'", "_");
            headerList.add(header);
        }
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            List<Object> values = new ArrayList<Object>();
            for (int i = 1; i < tokens.length; i++) {
                if (isBoolean)  {
                    if (tokens[i].equals("0"))
                        values.add(Boolean.FALSE);
                    else
                        values.add(Boolean.TRUE);
                }
                else
                    values.add(new Double(tokens[i]));
            }
            sampleToValues.put(tokens[0], values);
        }
        fu.close();
        //fileName = OVARIAN_DIR_NAME + "SampleToTOMNetworkModulesForAlteredGenes011110.arff";
        //fileName = OVARIAN_DIR_NAME + "SampleToTOMNetworkModulesForAlteredGenes011110_followup.arff";
        fileName = OVARIAN_DIR_NAME  + "SamplesToClustersForMutatedGenesSpectral_2Sample_040110.arff";
        fu.setOutput(fileName);
        printArffHeader(headerList, isBoolean);
        // Output data
        StringBuilder builder = new StringBuilder();
        for (String sample : sampleToPlatInterval.keySet()) {
            // A quick check to make sure this sample has expression value
            List<Object> values = sampleToValues.get(sample);
            if (values == null)
                continue;
            builder.setLength(0);
            Double platInt = sampleToPlatInterval.get(sample);
            if (platInt > 7.30)
                builder.append(false);
            else
                builder.append(true);
//            if (platInt <= 6.0) {
//                builder.append(true);
//            }
//            else if (platInt >= 9.0)
//                builder.append(false);
//            else
//                continue; // Don't want to output samples in between.
            for (Object value : values) {
                builder.append(",").append(value);
                // Some test code
                //                    if (Math.abs(1.0 - new Double(value.toString())) <= 1.0e-4)
                //                        builder.append(",true");
                //                    else
                //                        builder.append(",false");
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * This method is used to generate an arff for genes.
     * @throws Exception
     */
    @Test
    public void generateARFFFileForCorrGenes() throws Exception {
        PlatinumFreeIntervalScoreCalculator scoreCalculator = initScoreCalculator();
        Map<String, Map<String, Double>> geneExp = scoreCalculator.getGeneToExp();
        //List<String> genes = loadPlatCorrelatedGenes();
        // Gene HNRNPUL2 is not in the following list since it is not used in the validation data set
        String geneList = "PSG1, VNN1, TNFRSF10C, NOVA2, ZNF800, " +
        "NIPBL, BBC3, NUCB2, LIPN, ERGIC2, MLLT4, " +
        "FBXL20, CHAD, LOC641367, SLC4A8, FEN1, ZNF426, TARS, C12orf44";
        // The following list was generated from progress-free-time
        //String geneList = "DUSP23, BBC3, PSG1, TNFRSF10C, SLC4A8, EFNA1, ZFR, " +
        //      "VSTM2L, USP2, OR6M1, REG3A, CRELD1";
        String[] tmp = geneList.split(", ");
        List<String> genes = new ArrayList<String>(Arrays.asList(tmp));
        genes.retainAll(geneExp.keySet());
        System.out.println("Total genes: " + genes.size());
        Collections.sort(genes);
        List<Map<String, Double>> geneExpData = new ArrayList<Map<String,Double>>();
        List<String> attNames = new ArrayList<String>();
        // Try the first two graph components.
        for (int i = 0; i < genes.size(); i++) {
            String gene = genes.get(i);
            Map<String, Double> exp = scoreCalculator.getGeneToExp().get(gene);
            geneExpData.add(exp);
            attNames.add(gene);
            // Just a quick check
            //double score = scoreCalculator.calculateCorBetweenExpAndPlatInverval(exp, 
            //                                                                     scoreCalculator.getSampleToPlatinum());
            //System.out.println("Score for " + gene + ": " + score);
        }
        //String outputFileName = OVARIAN_DIR_NAME + "SampleExpToPlatIntervalGenesFroRbSurv20.arff";
        String outputFileName = OVARIAN_DIR_NAME + "DeleteMe.arff";
        generateARFF(scoreCalculator,
                     geneExpData,
                     attNames,
                     outputFileName);
    }
    
    private PlatinumFreeIntervalScoreCalculator initScoreCalculator() throws IOException {
        Map<String, Map<String, Double>> geneToExpData = loadGeneExp();
        Map<String, Double> sampleToPlatinum = loadSampleToPlatinumFreeInterval(true);
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        scoreCalculator.setGeneToExp(geneToExpData);
        scoreCalculator.setSampleToPlatinum(sampleToPlatinum);
        return scoreCalculator;
    }
    
    /**
     * This method is used to analyze correlation for FIs.
     * @throws Exception
     */
    @Test
    public void checkCorrelationOfGeneExpAndPlatForFIs() throws Exception {
        Map<String, Double> sampleToPlatinum = loadSampleToPlatinumFreeInterval();
        Map<String, Map<String, Double>> geneToExpData = loadGeneExp();
        Map<String, Double> geneToCor = calculateCorBetweenExpAndPlatInt(geneToExpData,
                                                                         sampleToPlatinum);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        System.out.println("Total FIs used: " + fis.size());
        int index = 0;
        final Map<String, Double> fiToCor = new HashMap<String, Double>();
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            Map<String, Double> exp1 = geneToExpData.get(gene1);
            Map<String, Double> exp2 = geneToExpData.get(gene2);
            if (exp1 == null || exp2 == null)
                continue;
            Map<String, Double> avgExp = averageExp(exp1, exp2);
            double cor = calculateCorBetweenExpAndPlatInverval(avgExp,
                                                               sampleToPlatinum);
            double cor1 = geneToCor.get(gene1);
            double cor2 = geneToCor.get(gene2);
            if ((cor > cor1 && cor > cor2) ||
                (cor < cor1 && cor < cor2)) {
                fiToCor.put(fi, cor);
            }
        }
        System.out.println("FIs with cor increased: " + fiToCor.size());
        List<String> fiList = new ArrayList<String>(fiToCor.keySet());
        Collections.sort(fiList, new Comparator<String>() {
            public int compare(String fi1, String fi2) {
                Double cor1 = fiToCor.get(fi1);
                Double cor2 = fiToCor.get(fi2);
                cor1 = Math.abs(cor1);
                cor2 = Math.abs(cor2);
                return cor2.compareTo(cor1);
            }
        });
        for (int i = 0; i < fiList.size(); i++) {
            String fi = fiList.get(i);
            Double cor = fiToCor.get(fi);
            System.out.println(fi + "\t" + cor);
            if (i > 100)
                break;
        }
    }
    
    private Map<String, Double> averageExp(Map<String, Double> exp1,
                                           Map<String, Double> exp2) {
        Map<String, Double> avg = new HashMap<String, Double>();
        for (String sample : exp1.keySet()) {
            Double value1 = exp1.get(sample);
            Double value2 = exp2.get(sample);
            if (value2 == null)
                continue;
            Double avgValue = (value1 + value2) / 2.0;
            avg.put(sample, avgValue);
        }
        return avg;
    }
    
    /**
     * This method is used to check correlation between gene expression and platinum interval.
     * @throws IOException
     */
    @Test
    public void checkCorrelationBetweenGeneExpAndPlatinum() throws IOException {
        Map<String, Double> sampleToPlatinum = loadSampleToPlatinumFreeInterval();
        Map<String, Map<String, Double>> geneToExpData = loadGeneExp();
        // Check correlation for each gene
        List<Double> platInterval = new ArrayList<Double>();
        List<Double> geneExp = new ArrayList<Double>();
        List<String> validSamples = new ArrayList<String>();
        final Map<String, Double> geneToCor = calculateCorBetweenExpAndPlatInt(geneToExpData, 
                                                                               sampleToPlatinum);
        List<String> geneList = new ArrayList<String>(geneToCor.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double cor1 = geneToCor.get(gene1);
                cor1 = Math.abs(cor1);
                Double cor2 = geneToCor.get(gene2);
                cor2 = Math.abs(cor2);
                return cor2.compareTo(cor1);
            }
        });
        // Print out the first 100 genes
        for (int i = 0; i < geneList.size(); i++) {
            String gene = geneList.get(i);
            Double cor = geneToCor.get(gene);
            if (Math.abs(cor) < 0.15)
                break;
            System.out.println(gene + "\t" + cor);
        }
        // Print out AMPH, which is included in the first network cluster
        Map<String, Double> amphGeneToExp = geneToExpData.get("AMPH");
        printGeneExpToPlatInterval(sampleToPlatinum,
                                   amphGeneToExp);
    }
    
    private List<String> loadPlatCorrelatedGenes() throws IOException {
        String fileName = OVARIAN_DIR_NAME + "CorrelationOfGeneExpAndPlatIntervalBasedOnMedian.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> genes = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double cor = new Double(tokens[1]);
            if (Math.abs(cor) >= 0.20)
                genes.add(tokens[0]);
        }
        fu.close();
        return genes;
    }
    
    /**
     * This method is used to annotate platinum correlated genes's pathways.
     * @param sample
     * @param sampleToValue
     * @return
     */
    @Test
    public void annotatePlatiumCorrelatedGenes() throws Exception {
        List<String> genes = loadPlatCorrelatedGenes();
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
//        annotator.annotateGenesWithFDR(genes);
//        System.out.println("\nAnnotate with one hop genes:");
//        annotator.setFDRThreshold(0.10);
//        annotator.annotateGenesWithFDROneHop(genes);
//        System.out.println("\nAnnotate with GO CC:");
//        annotator.annotateGenesUsingGOWithFDR(genes, "CC");
//        System.out.println("\nAnnotate with GO BP:");
//        annotator.annotateGenesUsingGOWithFDR(genes, "BP");
//        System.out.println("\nAnnotate with GO MF:");
//        annotator.annotateGenesUsingGOWithFDR(genes, "MF");
    }
    
    private double getExpressionValue(String sample, Map<String, Double> sampleToValue) {
        for (String expSample : sampleToValue.keySet()) {
            if (expSample.startsWith(sample)) {
                // Have to make sure the sample type is a cancer sample (sampleId < 10)
                String sub = expSample.substring(sample.length());
                int sampleId = Integer.parseInt(sub.substring(1, 3)); // get code 01 in -01A- 
                if (sampleId < 10)
                    return sampleToValue.get(expSample);
            }
        }
        return Double.NaN;
    }
    
    private boolean hasExpressionValue(String sample,
                                       Map<String, Double> sampleToExp) {
        for (String expSample : sampleToExp.keySet()) {
            if (expSample.startsWith(sample)) {
                // Have to make sure the sample type is a cancer sample (sampleId < 10)
                String sub = expSample.substring(sample.length());
                int sampleId = Integer.parseInt(sub.substring(1, 3)); // get code 01 in -01A- 
                if (sampleId < 10)
                    return true;
            }
        }
        return false;
    }
    
    /**
     * This method is used to load sample to platinum free interval.
     * @return
     * @throws IOException
     */
    public Map<String, Double> loadSampleToPlatinumFreeInterval() throws IOException {
        return loadSampleToPlatinumFreeInterval(false);
    }
    
    public Map<String, Double> loadSampleToSurvivalRate(boolean excludeLastFlowup) throws IOException {
        return new TCGAClinicalInformationLoader().loadSampleToSurvivalRate(excludeLastFlowup, CLINICAL_FILE);
    }
    
    public Map<String, Double> loadSampleToPlatinumFreeInterval(boolean excludeLastFlowUp) throws IOException {
        return new TCGAClinicalInformationLoader().loadSampleToPlatinumFreeInterval(excludeLastFlowUp, CLINICAL_FILE);
    }
    
    public Map<String, String> loadSampleToPrimaryTherapyOutcome() throws IOException {
        return new TCGAClinicalInformationLoader().loadSampleToPrimaryTherapyOutcome(CLINICAL_FILE);
    }
    
    @Test
    public void testLoadSampleToPlatinumFreeInterval() throws IOException {
        Map<String, Double> sampleToInterval = loadSampleToPlatinumFreeInterval();
        int index = 1;
        for (String sample : sampleToInterval.keySet()) {
            Double interval = sampleToInterval.get(sample);
            System.out.println(index + "\t" + sample + "\t" + interval);
            index ++;
        }
    }
    
    /**
     * Generate a FIs for genes.
     * @throws IOException
     */
    @Test
    public void buildNetworkForPlatCorrelatedGenes() throws IOException {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        List<String> genes = loadPlatCorrelatedGenes();
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        System.out.println("Total fis: " + fisInGenes.size());
        String fileName = OVARIAN_DIR_NAME + "FIsInPlatCorGenes.txt";
        fu.saveInteractions(fisInGenes, fileName);
    }
    
    /**
     * This method is used to do network module analysis for platinum correlated genes.
     * @throws Exception
     */
    @Test
    public void networkClusterGenesCorrelatedWithPlatinum() throws Exception {
        List<String> genes = loadPlatCorrelatedGenes();
        System.out.println("Total genes: " + genes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String fileName = OVARIAN_DIR_NAME + "PlatCorGenesNetworkClusters.txt";
        //List<Set<String>> clusters = clusterAnalyzer.cluster(genes);
        //clusterAnalyzer.outputNetworkClusters(clusters, fileName);
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(fileName);
        int index = 0;
        for (Set<String> cluster : clusters) {
            System.out.println("Cluster " + index + "\t" + cluster.size());
            index ++;
        }
        //clusterAnalyzer.annotateNetworkClusters(clusters, 5);
        // Check the correlation between the average expression values of network clusters
        // and samples.
        //Map<String, Map<String, Double>> geneToExpData = loadGeneExp();
        Map<String, Map<String, Double>> geneToExpData = loadTothillGeneExpData();
        //Map<String, Double> sampleToPlatInterval = loadSampleToPlatinumFreeInterval();
        Map<String, Double> sampleToPlatInterval = loadTothillSampleToRFI();
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        scoreCalculator.setGeneToExp(geneToExpData);
        scoreCalculator.setSampleToPlatinum(sampleToPlatInterval);
        scoreCalculator.setUseSampleIdDirectly(true);
        
        Map<String, Double> geneToCor = scoreCalculator.calculateScoreForGene();
        index = 0;
        for (Set<String> cluster : clusters) {
            Double score = scoreCalculator.calculateScore(cluster);
            System.out.println(index + "\t" + score + "\t" + cluster);
            for (String gene : cluster) {
                score = geneToCor.get(gene);
                System.out.println("\t" + gene + "\t" + score);
            }
            index ++;
        }
        // Want to print two lines
        // Print the first cluster
//        Map<String, Double> firstClusterToExp = clusterToExp.get(0);
//        printGeneExpToPlatInterval(sampleToPlatInterval,
//                                   firstClusterToExp);
    }

    private void printGeneExpToPlatInterval(Map<String, Double> sampleToPlatInterval,
                                            Map<String, Double> firstClusterToExp) {
        System.out.println("Sample\tExpression\tPlat_Interval");
        for (String sample : sampleToPlatInterval.keySet()) {
            if (hasExpressionValue(sample, firstClusterToExp)) {
                Double exp = getExpressionValue(sample, firstClusterToExp);
                Double interval = sampleToPlatInterval.get(sample);
                System.out.println(sample + "\t" + exp + "\t" + interval);
            }
        }
    }
    
    private Map<String, Double> calculateCorBetweenExpAndPlatInt(Map<?, Map<String, Double>> keyToExp,
                                                                 Map<String, Double> sampleToValue) {
        Map<String, Double> keyToCor = new HashMap<String, Double>();
        // Check correlation for each gene
        List<Double> platInterval = new ArrayList<Double>();
        List<Double> geneExp = new ArrayList<Double>();
        List<String> validSamples = new ArrayList<String>();
        for (Object key : keyToExp.keySet()) {
            Map<String, Double> sampleToExp = keyToExp.get(key);
            if (platInterval.size() == 0) {
                // Need to popup interval
                for (String sample : sampleToValue.keySet()) {
                    // Check if this sample has expression value
                    if (hasExpressionValue(sample, sampleToExp)) {
                        validSamples.add(sample);
                        platInterval.add(sampleToValue.get(sample));
                    }
                }
            }
            // Parse expression value
            geneExp.clear();
            for (String sample : validSamples) {
                double expValue = getExpressionValue(sample, sampleToExp);
                geneExp.add(expValue);
            }
            double cor = MathUtilities.calculatePearsonCorrelation(platInterval, geneExp);
            //System.out.println(key + "\t" + cor);
            keyToCor.put(key.toString(), cor);
        }
        return keyToCor;
    }
    
    private double calculateCorBetweenExpAndPlatInverval(Map<String, Double> sampleToExp,
                                                         Map<String, Double> sampleToInt) {
        return new PlatinumFreeIntervalScoreCalculator().calculateCorBetweenExpAndPlatInverval(sampleToExp, 
                                                                                               sampleToInt);
    }
    
    /**
     * This method is used to generate an ARFF for the Tothill data set.
     * @throws IOExecption
     */
    @Test
    public void generateARFFForTothillDatasets() throws IOException {
        Map<String, Map<String, Double>> geneExpData = loadTothillGeneExpData();
        Map<String, Double> sampleToRFI = loadTothillSampleToRFI();
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        scoreCalculator.setGeneToExp(geneExpData);
        scoreCalculator.setSampleToPlatinum(sampleToRFI);
        scoreCalculator.setUseSampleIdDirectly(true);
        //String fileName = TOTHILL_DIR + "SampleExpToRFIHieClustersGlobalZ.arff";
        //generateARFFForHiClusters(scoreCalculator, fileName);
        
        String fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatPos121809.txt";
        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
        List<GraphComponent> components = componentAnalyzer.loadComponentsInOrder(fileName);
        List<Map<String, Double>> componentExpData = new ArrayList<Map<String,Double>>();
        List<String> attNames = new ArrayList<String>();
        // Try the first two graph components.
        for (int i = 0; i < 1; i++) {
            GraphComponent comp = components.get(i);
            Set<String> genes = comp.getAllNodes();
            Map<String, Double> compExp = scoreCalculator.calculateAverage(genes);
            componentExpData.add(compExp);
            attNames.add("PosComponent" + i);
            // Just a quick check
            double score = scoreCalculator.calculateScore(comp);
            System.out.println("Score for component" + i + ": " + score);
        }
        // For negative part
        fileName = OVARIAN_DIR_NAME + "NetworkModulesFromSearchingForPlatNeg.txt";
        components = componentAnalyzer.loadComponentsInOrder(fileName);
        // Try the first two graph components.
        for (int i = 0; i < 1; i++) {
            GraphComponent comp = components.get(i);
            Set<String> genes = comp.getAllNodes();
            Map<String, Double> compExp = scoreCalculator.calculateAverage(genes);
            componentExpData.add(compExp);
            attNames.add("NegComponent" + i);
            // Just a quick check
            double score = scoreCalculator.calculateScore(comp);
            System.out.println("Score for component" + i + ": " + score);
        }
        if (true)
            return;
        String outputFileName = TOTHILL_DIR + "SampleExpToRFICompPosNegGlobalZ.arff";
        generateARFF(scoreCalculator,
                     componentExpData,
                     attNames,
                     outputFileName);
    }
    
    private Map<String, Double> loadTothillSampleToRFI() throws IOException {
        String fileName = TOTHILL_DIR + "SampleToRFIOrDeath.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Double> sampleToRFI = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0)
                break; // Get all data 
            //System.out.println(line);
            String[] tokens = line.split("\t");
            if (!isValidTothihllSample(tokens))
                continue;
            String sample = tokens[0];
            String value = tokens[12];
            sampleToRFI.put("X" + sample, 
                            new Double(value));
        }
        fu.close();
        return sampleToRFI;
    }
    
    protected Map<String, Double> loadTothillSampleToDeath() throws IOException {
        String fileName = TOTHILL_DIR + "SampleToRFIOrDeath.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Double> sampleToRFI = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0 || line.startsWith("#"))
                break; // Get all data 
            //System.out.println(line);
            String[] tokens = line.split("\t");
            if (!isValidTothihllSample(tokens))
                continue;
            String sample = tokens[0];
            String value = tokens[13];
            sampleToRFI.put("X" + sample, 
                            new Double(value));
        }
        fu.close();
        return sampleToRFI;
    }
    
    /**
     * Check if a tothill sample is a valid one 
     * @param tokens
     * @return
     */
    private boolean isValidTothihllSample(String[] tokens) {
        String sample = tokens[0];
        String value = tokens[12];
        if (value.length() == 0)
            return false;
        String primarySite = tokens[4];
        if (!primarySite.equals("OV"))
            return false; // Want to handle OV samples only
        // Make sure status value is not null
        String status = tokens[8];
        if (status.length() == 0)
            return false;
//        // Make sure platinum is used so that it can be compared with the TCGA data set
//        String pltx = tokens[9];
//        if (pltx.equals("N"))
//            return false;
        return true;
    }
    
    protected Map<String, String> loadTothillSampleToStatus() throws IOException {
        String fileName = TOTHILL_DIR + "SampleToRFIOrDeath.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, String> sampleToStatus = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0 || line.startsWith("#"))
                break; // Get all data 
            //System.out.println(line);
            String[] tokens = line.split("\t");
            if (!isValidTothihllSample(tokens))
                continue;
            String sample = tokens[0];
            String status = tokens[8];
            sampleToStatus.put("X" + sample, 
                               status);
        }
        fu.close();
        return sampleToStatus;
    }
    
    private Map<String, String> loadTothillGSMIdToSampleId() throws IOException {
        String fileName = TOTHILL_DIR + "GSMIDToSample.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, String> gsmId2SampleId = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            gsmId2SampleId.put(tokens[0], tokens[1]);
        }
        fu.close();
        return gsmId2SampleId;
    }
    
    @Test
    public void generateGSMTothillClinFile() throws IOException {
        Map<String, String> gsimToId = loadTothillGSMIdToSampleId();
        Map<String, String> idToGsim = new HashMap<String, String>();
        for (String gsim : gsimToId.keySet())
            idToGsim.put(gsimToId.get(gsim), gsim);
        String output = OVARIAN_DIR_NAME + "TothillClinInfo111811.txt";
        fu.setOutput(output);
        String input = OVARIAN_DIR_NAME + "TothillClinInfo102610.txt";
        fu.setInput(input);
        String line = fu.readLine();
        fu.printLine(line);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gsim = idToGsim.get(tokens[0]);
            fu.printLine(gsim + "\t" + tokens[1] + "\t" + tokens[2]);
        }
        fu.close();
    }
    
    @Test
    public void generateTothillSampleToGeneExpScoreForModule4() throws Exception {
        // Need to load gene sets
        String clusterFileName = OVARIAN_DIR_NAME + "ExomeMutationModules_3More_091310.txt";
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> modules = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        // Get want to get the fifth module
        Set<String> module = modules.get(4);
        System.out.println("Genes in module: " + module.size() + ": " + module);
        Map<String, Map<String, Double>> geneToSampleToExp = loadTothillGeneExpData();
        //Map<String, Map<String, Double>> geneToSampleToExp = new GSEADataHelper().loadGeneExp(TOTHILL_DIR + "Tothill_TCGA_plus2_sample_gene_z.txt");
        module.retainAll(geneToSampleToExp.keySet());
        System.out.println("Shared genes: " + module.size());
        Map<String, String> sampleToStatus = loadTothillSampleToStatus();
        Map<String, Double> sampleToSurvival = loadTothillSampleToDeath();
        Map<String, List<Double>> sampleToValues = new HashMap<String, List<Double>>();
        for (String gene : module) {
            Map<String, Double> sampleToValue = geneToSampleToExp.get(gene);
            if (sampleToValue == null)
                continue;
            for (String sample : sampleToValue.keySet()) {
                Double value = sampleToValue.get(sample);
                List<Double> values = sampleToValues.get(sample);
                if (values == null) {
                    values = new ArrayList<Double>();
                    sampleToValues.put(sample, values);
                }
                values.add(value);
            }
        }
        System.out.println("Sample\tType\tAverage\tStatus\tSurvival");
        double cutoff = 3.685d;
        for (String sample : sampleToStatus.keySet()) {
            List<Double> values = sampleToValues.get(sample);
            if (values == null)
                continue;
            Double average = MathUtilities.calculateMean(values);
            String status = sampleToStatus.get(sample);
            Double survival = sampleToSurvival.get(sample);
            if (average > cutoff)
                System.out.println(sample + "\t1\t" + average + "\t" + status + "\t" + survival);
            else if (average < cutoff)
                System.out.println(sample + "\t-1\t" + average + "\t" + status + "\t" + survival);
        }
    }
    
    protected Map<String, Map<String, Double>> loadTothillGeneExpData() throws IOException {
        Map<String, String> gsmToSample = loadTothillGSMIdToSampleId();
        Map<String, Map<String, Double>> geneToData = new HashMap<String, Map<String,Double>>();
//        String fileName = TOTHILL_DIR + "Tothill_TCGA_plus2_z.txt";
//        String fileName = TOTHILL_DIR + "Tothill_TCGA_plus2_sample_gene_z.txt";
        //String fileName = TOTHILL_DIR + "Tothill_TCGA_plus2_global_median_centered.txt";
//        String fileName = TOTHILL_DIR + "GSE9899_z_Gene_Exp_111811.txt";
//        String fileName = TOTHILL_DIR + "GSE9899_Gene_Exp_111811.txt";
//        String fileName = TOTHILL_DIR + "GSE9899_GSE9891_Gene_Exp_111811.txt";
        String fileName = TOTHILL_DIR + "GSE9891_z_Gene_Exp_111811.txt";
        fu.setInput(fileName);
        // Sample list
        String line = fu.readLine();
        List<String> gsmIdList = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++)
            gsmIdList.add(tokens[i]);
        // Need to convert gsmId to sampleId
        List<String> sampleList = new ArrayList<String>();
        for (String gsmId : gsmIdList) {
            String sampleId = gsmToSample.get(gsmId);
            sampleList.add(sampleId);
        }
//        List<String> sampleList = new ArrayList<String>(gsmIdList);
        // Starting parsing
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // The first one is gene name
            String gene = tokens[0];
            Map<String, Double> sampleToValue = new HashMap<String, Double>();
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA"))
                    continue; // Just escape these values
                String sample = sampleList.get(i - 1);
                Double value = new Double(tokens[i]);
                sampleToValue.put(sample, value);
            }
            geneToData.put(gene, sampleToValue);
        }
        fu.close();
        return geneToData;
    }
    
    /**
     * This method is used to generate an gene expression data set plus recurrence-free information
     * based on genes identified from rbsurv R package.
     * @throws Exception
     */
    @Test
    public void generateTothillDataGeneExpDataForRbSurvGenes() throws Exception {
        Map<String, Map<String, Double>> geneToExp = loadTothillGeneExpData();
        // The following list was generated based on platinum-free-interval
        String geneList = "PSG1, VNN1, TNFRSF10C, NOVA2, ZNF800, " +
                          "NIPBL, BBC3, HNRNPUL2, NUCB2, LIPN, ERGIC2, MLLT4, " +
                          "FBXL20, CHAD, LOC641367, SLC4A8, FEN1, ZNF426, TARS, C12orf44";
        // The following list was generated from progress-free-time
//        String geneList = "DUSP23, BBC3, PSG1, TNFRSF10C, SLC4A8, EFNA1, ZFR, " +
//        		          "VSTM2L, USP2, OR6M1, REG3A, CRELD1";
        String[] tmp = geneList.split(", ");
        List<String> selectedGenes = Arrays.asList(tmp);
        for (String gene : selectedGenes)
            System.out.print("\"" + gene + "\", ");
        System.out.println();
        if (true)
            return;
        // Check how many genes are shared between these two data sets
        Set<String> shared = InteractionUtilities.getShared(selectedGenes, geneToExp.keySet());
        System.out.println("Shared genes: " + shared.size());
        // Check genes not in the tothill
        Set<String> copy = new HashSet<String>(selectedGenes);
        copy.removeAll(shared);
        System.out.println("Genes not in the tothill data set: " + copy);
        // The following are used to handle the clin data 
        Map<String, String> sampleToStatus = loadTothillSampleToStatus();
        Map<String, Double> sampleToRFI = loadTothillSampleToRFI();
        String outFileName = OVARIAN_DIR_NAME + "TothillGeneExpForRbsurvGenesForProgreeFree.txt";
        fu.setOutput(outFileName);
        // Generate a header
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        List<String> sharedList = new ArrayList<String>(shared);
        Collections.sort(sharedList);
        for (String gene : sharedList)
            builder.append("\t").append(gene);
        builder.append("\tSTATUS");
        builder.append("\tRFI");
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String sample : sampleToRFI.keySet()) {
            String status = sampleToStatus.get(sample);
            Double rfi = sampleToRFI.get(sample);
            builder.append(sample);
            for (String gene : sharedList) {
                Map<String, Double> sampleToExp = geneToExp.get(gene);
                Double exp = sampleToExp.get(sample);
                builder.append("\t").append(exp);
            }
            if (status.equals("PF"))
                builder.append("\t0");
            else
                builder.append("\t1");
            builder.append("\t").append(rfi);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    
    @Test
    public void ztransformExpDataBasedOnGene() throws Exception {
        Map<String, Map<String, Double>> geneToSampleToExp = loadTothillGeneExpData();
        DescriptiveStatistics stat = new DescriptiveStatistics();
        StringBuilder builder = new StringBuilder();
        String outFile = TOTHILL_DIR + "Tothill_TCGA_plus2_sample_gene_z.txt";
        fu.setOutput(outFile);
        boolean isHeaderDone = false;
        List<String> sampleList = new ArrayList<String>();
        for (String gene : geneToSampleToExp.keySet()) {
            Map<String, Double> sampleToExp = geneToSampleToExp.get(gene);
            if (!isHeaderDone) {
                isHeaderDone = true;
                // Generate a header
                sampleList.addAll(sampleToExp.keySet());
                Collections.sort(sampleList);
                builder.append("Gene");
                for (String sample : sampleList)
                    builder.append("\t").append(sample);
                fu.printLine(builder.toString());
                builder.setLength(0);
            }
            if (isHeaderDone) {
                // Make sure all have the same length
                if (sampleToExp.size() != sampleList.size())
                    throw new IllegalStateException("Sample has different length for gene " + gene);
            }
            for (String sample : sampleList) {
                Double value = sampleToExp.get(sample);
                stat.addValue(value);
            }
            double mean = stat.getMean();
            double sd = stat.getStandardDeviation();
            // Did a transformation
            builder.append(gene);
            for (String sample : sampleList) {
                Double value = sampleToExp.get(sample);
                value = (value - mean) / sd;
                builder.append("\t").append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
            stat.clear();
        }
        fu.close();
    }
    
    /**
     * This helper method is used to divide samples into groups with similar survival rate distribution.
     * @param sampleToGenes
     * @param sampleToClinInfo
     * @param fold
     * @return
     * @throws Exception
     */
    private List<Set<String>> stratifiedCrossValidationSplit(Map<String, Set<String>> sampleToGenes,
                                                            Map<String, TCGAClinicalInfo> sampleToClinInfo,
                                                            int fold,
                                                            int binNumber) throws Exception {
        // Used for discretizing
        final Map<String, Double> sampleToRate = new HashMap<String, Double>();
        for (String sample : sampleToGenes.keySet()) {
            TCGAClinicalInfo info = sampleToClinInfo.get(sample);
            if (info == null || info.getOsEvent() == null || info.getOsDuration() == null)
                continue;
            sampleToRate.put(sample, info.getOsDuration());
        }
        List<String> sampleList = new ArrayList<String>(sampleToRate.keySet());
        Collections.sort(sampleList, new Comparator<String>() {
            public int compare(String sample1, String sample2) {
                Double value1 = sampleToRate.get(sample1);
                Double value2 = sampleToRate.get(sample2);
                return value2.compareTo(value1);
            }
        });
        // Create bins based on rate distribution
        List<Set<String>> sampleBins = new ArrayList<Set<String>>();
        int numberInBin = sampleList.size() / binNumber;
        for (int i = 0; i < binNumber; i++) {
            // Pick up the top numberInBin
            Set<String> sampleBin = new HashSet<String>(sampleList.subList(0, numberInBin));
            sampleBins.add(sampleBin);
            sampleList.removeAll(sampleBin);
        }
        // Check how many samples for each group
        List<Set<String>> dividedSamples = new ArrayList<Set<String>>();
        for (int i = 0; i < fold; i++) {
            dividedSamples.add(new HashSet<String>());
        }
        RandomData randomizer = new RandomDataImpl();
        for (Set<String> sampleBin : sampleBins) {
            // Divide samples in each bin into n fold
            int size = sampleBin.size() / fold;
            for (int i = 0; i < fold; i++) {
                Set<String> randomSample = MathUtilities.randomSampling(sampleBin, 
                                                                        size, 
                                                                        randomizer);
                sampleBin.removeAll(randomSample);
                Set<String> group = dividedSamples.get(i);
                group.addAll(randomSample);
            }
        }
        // The following code is used to check the distribution of samples
        int index = 0;
        for (Set<String> group : dividedSamples) {
            SummaryStatistics stat = new SummaryStatistics();
            for (String sample : group) {
                TCGAClinicalInfo info = sampleToClinInfo.get(sample);
                stat.addValue(info.getOsDuration());
            }
            System.out.println("Group " + index + ": " + stat.getMean() + ", " +
                               stat.getStandardDeviation());
            index ++;
        }
        return dividedSamples;
    }
    
    /**
     * This method is used to do cross-validation test for network modules that are related to clinical information.
     * @throws Exception
     */
    @Test
    public void crossValidationTestOnClinInfoForModules() throws Exception {
        // Load clinical information
        Map<String, TCGAClinicalInfo> sampleToClinInfo = new TCGAClinicalInformationLoader().loadClinicalInfo(CLINICAL_FILE);
        
        //Map<String, Set<String>> sampleToGenes = loadSampleToNonSynonymousMutatedGenes(false);
        Map<String, Set<String>> sampleToGenes = loadSampleToAlteredGenes();
        Set<String> upperSamples = new HashSet<String>();
        Set<String> lowerSamples = new HashSet<String>();
        divideSamplesIntoUpperLowerOnSurvival(upperSamples,
                                              lowerSamples,
                                              sampleToGenes);
        // Get the divided samples
        int fold = 5;
//        List<Set<String>> dividedSamples = randomCrossValidationSplit(sampleToGenes, fold);
        List<Set<String>> dividedSamples = stratifiedCrossValidationSplit(sampleToGenes, 
                                                                          sampleToClinInfo, 
                                                                          fold, 
                                                                          20);
        // Check the number of samples
        for (Set<String> randomSample : dividedSamples)
            System.out.println("Random sample: " + randomSample.size());
        
        SpectralPartitionNetworkCluster clustering = new SpectralPartitionNetworkCluster();
        List<List<Set<String>>> dividedClusters = new ArrayList<List<Set<String>>>();
        for (int i = 0; i < fold; i++) {
            System.out.println("\nValidation sample: " + i);
            Set<String> validationSamples = dividedSamples.get(i);
            Map<String, Set<String>> copy = new HashMap<String, Set<String>>(sampleToGenes);
            copy.keySet().removeAll(validationSamples);
            // Need to make a copy to avoid any changes to the original set
            for (String sample : copy.keySet()) {
                Set<String> set = copy.get(sample);
                copy.put(sample, new HashSet<String>(set));
            }
            filterSampleToGenes(copy, 2);
            Set<String> genes = InteractionUtilities.grepAllGenes(copy);
            System.out.println("Total genes: " + genes.size());
            List<Set<String>> clusters = clustering.clusterGenes(genes);
            dividedClusters.add(clusters);
            System.out.println("Validation samples:");
            fisherExactTest(upperSamples, 
                            lowerSamples, 
                            validationSamples,
                            sampleToGenes, 
                            clusters);
            System.out.println("Training samples:");
            fisherExactTest(upperSamples,
                            lowerSamples,
                            copy.keySet(),
                            sampleToGenes,
                            clusters);
            // Generate a sample to module file
//            String fileName = OVARIAN_DIR_NAME + "CrossValidationSampleAlteredGenesToModules_" + fold + "fold_" + i + ".txt";
            String fileName = OVARIAN_DIR_NAME + "StratifiedCrossValidation2SampleAlteredGenesToModules_" + fold + "fold_" + i + ".txt";
            outputSampleToModules(clusters, 
                                  sampleToGenes,
                                  validationSamples,
                                  fileName);
            // There is a bug in javastat.jar file. Java log rank actually cannot work.
//            // The following code is used to log-rank analaysis
//            for (int j = 0; j < clusters.size(); j++) {
//                Set<String> cluster = clusters.get(j);
//                if (cluster.size() < 10)
//                    continue;
//                System.out.println("Module " + j + ": " + cluster.size());
//                System.out.println("Validation results:");
//                logRankAnalysis(sampleToClinInfo, 
//                                cluster, 
//                                sampleToGenes,
//                                validationSamples);
//                System.out.println("Training results:");
//                logRankAnalysis(sampleToClinInfo, 
//                                cluster, 
//                                sampleToGenes,
//                                copy.keySet());
//            }
        }
//        System.out.println("\nClusters overlapping analysis:");
//        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
//        filterSampleToGenes(sampleToGenes, 2);
//        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
//        for (int i = 0; i < dividedClusters.size() - 1; i++) {
//            List<Set<String>> clusters1 = dividedClusters.get(i);
//            for (int j = i + 1; j < dividedClusters.size(); j++) {
//                List<Set<String>> clusters2 = dividedClusters.get(j);
//                System.out.println("\n" + i + " vs " + j);
//                clusterAnalyzer.checkNetworkModuleOverlapping(clusters1, clusters2, 4, 0.01, allGenes);
//            }
//        }
    }
    
    private void outputSampleToModules(List<Set<String>> clusters,
                                       Map<String, Set<String>> sampleToGenes,
                                       Set<String> validationSamples,
                                       String fileName) throws IOException {
        int size = 4;
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample\tValidation");
        for (int i = 0; i < clusters.size(); i++) {
            if (clusters.get(i).size() < 4)
                break;
            builder.append("\tModule").append(i);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            builder.append(sample).append("\t");
            if (validationSamples.contains(sample))
                builder.append(1);
            else
                builder.append(0);
            for (int i = 0; i < clusters.size(); i++) {
                Set<String> cluster = clusters.get(i);
                if (cluster.size() < size)
                    break;
                builder.append("\t");
                if (InteractionUtilities.isShared(cluster, genes))
                    builder.append(1);
                else
                    builder.append(0);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private void fisherExactTest(Set<String> upperSamples,
                                 Set<String> lowerSamples,
                                 Set<String> validationSamples,
                                 Map<String, Set<String>> sampleToGenes,
                                 List<Set<String>> clusters) {
        Set<String> longSamples = InteractionUtilities.getShared(upperSamples, validationSamples);
        Set<String> shortSamples = InteractionUtilities.getShared(lowerSamples, validationSamples);
        System.out.println("Total upper samples: " + longSamples.size());
        System.out.println("Total lower samples: " + shortSamples.size());
        // For Fisher exact test
        System.out.println("Module\tSize\tSamplesInUpper\tSamplesInLower\tTwo_tail_P\tLeft_P\tRight_P");
        int a, b;
        FisherExact fisherTest = new FisherExact(sampleToGenes.size());
        for (int j = 0; j < clusters.size(); j++) {
            Set<String> cluster = clusters.get(j);
            a = b = 0;
            for (String sample : longSamples) {
                Set<String> genes1 = sampleToGenes.get(sample);
                Set<String> shared = InteractionUtilities.getShared(genes1, cluster);
                if (shared.size() > 0)
                    a ++;
            }
            for (String sample : shortSamples) {
                Set<String> genes1 = sampleToGenes.get(sample);
                Set<String> shared = InteractionUtilities.getShared(genes1, cluster);
                if (shared.size() > 0)
                    b ++;
            }
            double pvalue = fisherTest.getTwoTailedP(a, b, longSamples.size() - a, shortSamples.size() - b); 
            double left = fisherTest.getLeftTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            double right = fisherTest.getRightTailedP(a, b, longSamples.size() - a, shortSamples.size() - b);
            
            System.out.println("Module " + j + "\t" + cluster.size() + "\t" + a + "\t" + b + "\t" + String.format("%.3f", pvalue) + "\t" + 
                               String.format("%.3f", left) + "\t" + String.format("%.3f", right));
        }
    }
    
    private void logRankAnalysis(Map<String, TCGAClinicalInfo> sampleToClinInfo,
                                 Set<String> module,
                                 Map<String, Set<String>> sampleToGenes,
                                 Set<String> samples) {
        // Split samples into two groups: one for samples having mutated genes, another not
        Set<String> samplesInModules = new HashSet<String>();
        Set<String> samplesNotInModule = new HashSet<String>();
        for (String sample : samples) {
            Set<String> genes = sampleToGenes.get(sample);
            if (InteractionUtilities.isShared(genes, module))
                samplesInModules.add(sample);
            else
                samplesNotInModule.add(sample);
        }
        List<Double> time1 = new ArrayList<Double>();
        List<Double> censor1 = new ArrayList<Double>();
        for (String sample : samplesInModules) {
            TCGAClinicalInfo info = sampleToClinInfo.get(sample);
            if (info == null || info.getOsDuration() == null || info.getOsEvent() == null)
                continue;
            time1.add(info.getOsDuration());
            if (info.getOsEvent())
                censor1.add(1.0d);
            else
                censor1.add(0.0d);
        }
        List<Double> time2 = new ArrayList<Double>();
        List<Double> censor2 = new ArrayList<Double>();
        for (String sample : samplesNotInModule) {
            TCGAClinicalInfo info = sampleToClinInfo.get(sample);
            if (info == null|| info.getOsDuration() == null || info.getOsEvent() == null)
                continue;
            time2.add(info.getOsDuration());
            if (info.getOsEvent())
                censor2.add(1.0d);
            else
                censor2.add(0.0d);
        }
        System.out.println(""+ time1 + censor1 + time2 + censor2);
//        double pvalue = MathUtilities.logRankSurvivalTest(time1, censor1, time2, censor2);
//        System.out.println("pvalue from log-rank: " + pvalue);
    }
    
    /**
     * This method is used to check overlapping of correlated genes between TCGA data set
     * and Tothill data set.
     * @throws Exception
     */
    @Test
    public void analyzeCorGeneOverlapping() throws Exception {
        Map<String, Map<String, Double>> geneToExpData = loadTothillGeneExpData();
        Map<String, Double> sampleToPlatInterval = loadTothillSampleToRFI();
        PlatinumFreeIntervalScoreCalculator scoreCalculator = new PlatinumFreeIntervalScoreCalculator();
        scoreCalculator.setGeneToExp(geneToExpData);
        scoreCalculator.setSampleToPlatinum(sampleToPlatInterval);
        scoreCalculator.setUseSampleIdDirectly(true);
        
        Map<String, Double> geneToCor = scoreCalculator.calculateScoreForGene();
        // Get a list of genes with cor >= 0.15
        List<String> tothillGenes = new ArrayList<String>();
        for (String gene : geneToCor.keySet()) {
            Double cor = geneToCor.get(gene);
            if (Math.abs(cor) >= 0.20)
                tothillGenes.add(gene);
        }
        // Compare with tcga genes
        List<String> tcgaGenes = loadPlatCorrelatedGenes();
        // Get the total genes
        Map<String, Map<String, Double>> tcgaGeneExp = loadGeneExp();
        Set<String> sharedTotal = new HashSet<String>(tcgaGeneExp.keySet());
        sharedTotal.retainAll(geneToCor.keySet());
        
        Set<String> shared = new HashSet<String>(tcgaGenes);
        shared.retainAll(tothillGenes);
        
        System.out.println("Total shared genes: " + sharedTotal.size());
        System.out.println("Shared: " + shared.size());
        System.out.println("Total tothill genes: " + tothillGenes.size());
        System.out.println("Total TCGA genes: " + tcgaGenes.size());
        
        double pvalue = MathUtilities.calculateHypergeometricPValue(sharedTotal.size(), 
                                                                    tcgaGenes.size(),
                                                                    tothillGenes.size(), 
                                                                    shared.size());
        System.out.println("pvalue from hypergeometric: " + pvalue);
    }
    
}
