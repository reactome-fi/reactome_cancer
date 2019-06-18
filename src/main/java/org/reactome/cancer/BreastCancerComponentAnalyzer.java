/*
 * Created on Jun 27, 2008
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;
import org.gk.util.StringUtils;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.r3.EntrezGeneAnalyzer;
import org.reactome.r3.graph.GraphComponent;
import org.reactome.r3.graph.GraphComponentSearchEngine;
import org.reactome.r3.graph.MCLResultsAnalyzer;
import org.reactome.r3.graph.MutualInformationScoreCalculator;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

import cern.colt.list.ObjectArrayList;
 /** The major functions of
 * this class are refactored from class BreastCancerArrayDataAnalyser.
 * @author wgm
 *
 */
public class BreastCancerComponentAnalyzer {
    private final String dirName = "datasets/BreastCancer/";
    private FileUtility fu = new FileUtility();
    private Set<String> totalGenes = null;
    private Set<String> allFIs = null;
    private final double MAX_P_VALUE = 1.0E-5;
    private BreastCancerArrayDataAnalyzer dataAnalyzer = new BreastCancerArrayDataAnalyzer();
    private final int CUT_OFF_OF_COMP_SIZE = 4;
    
    public BreastCancerComponentAnalyzer() {
    }
    
    @Test
    public void checkComponentsInCompleteGraph() throws Exception {
        MutualInformationScoreCalculator calculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        calculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        calculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        // Want to use the biggest component only
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        List<GraphComponent> components = loadGSE2034Comonents();
        for (int i = 0; i < components.size(); i++) {
            GraphComponent comp = components.get(i);
            Set<String> nodes = comp.getAllNodes();
            // Just try one
            for (String seed : nodes) {
                GraphComponent newComp = searchComponentInCompleteGraph(seed, ids, calculator);
                System.out.println("Old Score: " + comp.getScore());
                System.out.println("New Score: " + newComp.getScore());
            }
            break;
        }
    }
    
    private GraphComponent searchComponentInCompleteGraph(String seed,
                                                          Set<String> ids,
                                                          MutualInformationScoreCalculator miCalculator) throws Exception {
        Set<String> copy = new HashSet<String>(ids);
        copy.remove(seed);
        GraphComponent comp = new GraphComponent();
        comp.addNode(seed);
        comp.setScoreCalculator(miCalculator);
        while (copy.size() > 0) {
            double oldScore = miCalculator.calculateScore(comp);
            int oldSize = comp.getAllNodes().size();
            for (String id : copy) {
                comp.addNode(id);
                double newScore = miCalculator.calculateScore(comp);
                if (newScore > oldScore) {
                    //System.out.println(comp.getAllNodes().size() + ": " + newScore);
                    oldScore = newScore;
                }
                else 
                    comp.removeNode(id);
            }
            if (comp.getAllNodes().size() == oldSize)
                break;
            copy.removeAll(comp.getAllNodes());
            break;
        }
        return comp;
    }
    
    @Test
    public void checkPermutatedComps() throws Exception {
        String compFileName = dirName + "Components2_1_GSE2023_Sample_Permutated.txt";
        List<GraphComponent> components = loadComponentsInOrder(compFileName);
        MutualInformationScoreCalculator calculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        calculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        calculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        System.out.println("Comp\tOldValue\tNewValue");
        for (GraphComponent comp : components) {
            Double oldValue = comp.getScore();
            Double newValue = calculator.calculateScore(comp);
            System.out.println(comp.getId() + "\t" + oldValue + "\t" + newValue);
        }
    }
    
    /**
     * This method is used to generate an arff file for weka to do classifying.
     * @throws Exception
     */
    @Test
    public void generateArffFileForComponents() throws Exception {
        //List<GraphComponent> filteredComps = loadClusters();
        //List<GraphComponent> filteredComps = loadMCLClusters();
        //filterRedundantComponents(filteredComps);
        //System.out.println("After redundant removed: " + filteredComps.size());
        // Need to get the data
        // Need a MutualInformationScoreCalculator
        MutualInformationScoreCalculator calculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        // Wang et al
        calculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        calculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
//        //String fileName = dirName + "Components2_1_GSE2023_Sample_Permutated_cutoff_10.arff";
        String fileName = dirName + "FIsWithWeightFromTTestMCLClustersGSE2034_I80_091509.arff";
//        // GSE4922
//        calculator.setSampleToPhenotypeInfo(loadGSE4922PatientInfo());
//        calculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore.txt");
//        String fileName = dirName + "Components2_1_GSE4922_NO_P_I_cutoff_10_GSE4922.arff";        
//        // For van de Vijver data set
        //calculator.setSampleToPhenotypeInfo(dataAnalyzer.loadNejmPatientInfo());
        //calculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
        //String fileName = dirName + "Components2_1_Nejm_NO_P_I_cutoff_11_Nejm.arff";
        //String fileName = dirName + "Components2_1_GSE2034_062609_Nejm.arff";
        String clusterFileName = dirName + "FIsWithWeightFromTTestMCLClustersGSE2034_I80_091509.txt";
        List<GraphComponent> filteredComps = loadMCLClusterFromTScore(clusterFileName, 
                                                                 calculator,
                                                                 2.0);
        List<String> samples = calculator.getSamples();
        Map<String, String> sampleToType = calculator.getSampleToType();
        Map<GraphComponent, List<Double>> compToValues = new HashMap<GraphComponent, List<Double>>();
        for (GraphComponent comp : filteredComps) {
            List<Double> values = calculator.averageValues(comp);
            compToValues.put(comp, values);
        }
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("@relation BreastCancerType\n\n");
        builder.append("@attribute CellType {true, false}\n");
        for (GraphComponent comp : filteredComps) {
            builder.append("@attribute Comp").append(comp.getId());
            builder.append(" NUMERIC\n");
        }
        builder.append("\n@data");
        fu.printLine(builder.toString());
        builder.setLength(0);
        // Output the table
        for (int i = 0; i < samples.size(); i++) {
            String sample = samples.get(i);
            String type = sampleToType.get(sample);
            if (type.equals("0"))
                type = "false";
            else
                type = "true";
            builder.append(type);
            for (GraphComponent comp : filteredComps) {
                List<Double> values = compToValues.get(comp);
                builder.append(",").append(values.get(i));
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private List<GraphComponent> loadMCLClusters() throws IOException {
        List<GraphComponent> clusters = new ArrayList<GraphComponent>();
        String fileName = dirName + "AverageFIWeightsForMCLClusters.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double score = Double.parseDouble(tokens[1]);
            if (score < 0.23)
                break;
            String[] genes = tokens[5].split(",");
            GraphComponent comp = new GraphComponent();
            comp.setId(Integer.parseInt(tokens[0]));
            for (String gene : genes)
                comp.addNode(gene);
            clusters.add(comp);
        }
        fu.close();
        return clusters;
    }

    private List<GraphComponent> loadClusters() throws IOException {
        double miCutOff = 0.10;
        //double miCutOff = 0.11;
        //double miCutOff = 0.09;
        //double miCutOff = 0.125d;
        //double miCutOff = 0.136d;
        //String compFileName = dirName + "Components2_1_GSE2023_Sample_Permutated.txt";
        String compFileName = dirName + "Components2_1_GSE2034_062609.txt";
        List<GraphComponent> components = loadComponentsInOrder(compFileName);
        List<GraphComponent> filteredComps = filterComponents(miCutOff,
                                                              components);
        System.out.println("Total components: " + components.size());
        System.out.println("Filtered components: " + filteredComps.size());
        return filteredComps;
    }
    
    @Test
    public void analyzeGraphCompFIDataset() throws Exception {
        List<GraphComponent> gse2034Comps = loadGSE2034Comonents();
        List<GraphComponent> gse4922Comps = loadGSE4922Components();
        Map<String, Double> pairToPValue = calculatePValuesForTwoComponentSets(gse2034Comps,
                                                                               gse4922Comps);
        String fileName = dirName + "GraphComponentFIComp.txt";
        //String fileName = dirName + "GraphComponentFICompRandom.txt";
        fu.setInput(fileName);
        Map<String, Set<String>> genesInGse2034 = new HashMap<String, Set<String>>();
        Map<String, Set<String>> genesInGse4922 = new HashMap<String, Set<String>>();
        String line = fu.readLine();
        boolean isGse2034 = false;
        boolean isGse4922 = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#GSE2034")) {
                isGse2034 = true;
                isGse4922 = false;
                continue;
            }
            else if (line.startsWith("#GSE4922")) {
                isGse4922 = true;
                isGse2034 = false;
                continue;
            }
            else if (line.startsWith("#NEJM")) {
                isGse4922 = false;
                isGse2034 = false;
                continue;
            }
            String[] tokens = line.split("\t");
            double pvalue = Double.parseDouble(tokens[3]);
            if (pvalue > 0.0005)
                continue;
            String comp = tokens[0];
            int index = comp.indexOf(" ");
            comp = comp.substring(index + 1);
            String geneNames = tokens[4];
            geneNames = geneNames.substring(1, geneNames.length() - 1);
            tokens = geneNames.split(",");
            Set<String> genes = new HashSet<String>();
            for (String token : tokens)
                genes.add(token.trim());
            if (isGse2034)
                genesInGse2034.put(comp, genes);
            else if (isGse4922)
                genesInGse4922.put(comp, genes);
        }
        fu.close();
        // Check pair-wise p-values
        int total = 45;
        Set<String> filteredGse2034 = new HashSet<String>();
        Set<String> filteredGse4922 = new HashSet<String>();
        for (String comp1: genesInGse2034.keySet()) {
            Set<String> set1 = genesInGse2034.get(comp1);
            for (String comp2 : genesInGse4922.keySet()) {
                String key = comp1 + " " + comp2;
                if (!pairToPValue.containsKey(key))
                    continue;
                Set<String> set2 = genesInGse4922.get(comp2);
                int shared = sharedElements(set1, set2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(total, 
                                                                            set1.size(), 
                                                                            set2.size(), 
                                                                            shared);
                if (pvalue > MAX_P_VALUE)
                    continue;
                Double compPValue = pairToPValue.get(key);
                System.out.println(key + "\t" + pvalue + "\t" + compPValue); 
                filteredGse2034.add(comp1);
                filteredGse4922.add(comp2);
            }
        }
        System.out.println("Filtered components in GSE2034: " + 
                           filteredGse2034.size() + ": " + 
                           filteredGse2034);
        System.out.println("Filtered components in GSE4922: " + 
                           filteredGse4922.size() + ": " +
                           filteredGse4922);
    }
    
    @Test
    public void calculatePValueForComponentsInDatasets() throws Exception {
        // From Wang et al
        //String compFileName = dirName + "Components2_1_Norm.txt";
        //List<GraphComponent> gse2034Comps = loadComponentsInOrder(compFileName);
        //
        // For GSE4922
        String compFileName = dirName + "Components2_1_GSE4922_Norm.txt";
        List<GraphComponent> gse4922Comps = loadComponentsInOrder(compFileName);
        List<GraphComponent> comps = gse4922Comps;
        int lastIndex = 51;
        Map<String, Double> pairToPValue = new HashMap<String, Double>();
        for (int i = 0; i < lastIndex - 1; i++) {
            GraphComponent comp1 = comps.get(i);
            Set<String> ids1 = comp1.getAllNodes();
            for (int j = i + 1; j < lastIndex; j++) {
                GraphComponent comp2 = comps.get(j);
                Set<String> ids2 = comp2.getAllNodes();
                int shared = sharedElements(ids1, ids2);
                double p = MathUtilities.calculateHypergeometricPValue(9500, 
                                                                       ids1.size(), 
                                                                       ids2.size(), 
                                                                       shared);
                if (p > MAX_P_VALUE)
                    continue;
                String key = i + " " + j;
                System.out.println(key +": " + p);
                pairToPValue.put(key, p);
            }
        }
    }
    
    @Test
    public void calculatePValuesAmongComponents() throws Exception {
        List<GraphComponent> gse2034Comps = loadGSE2034Comonents();
        List<GraphComponent> gse4922Comps = loadGSE4922Components();
        List<GraphComponent> nejmComps = loadNejmComponents();
        System.out.println("Components between GSE2034 to nejm:");
        Map<String, Double> pairToValue = calculatePValuesForTwoComponentSets(gse2034Comps, 
                                                                              nejmComps);
        List<String> keyList = new ArrayList<String>(pairToValue.keySet());
        Collections.sort(keyList);
        for (String key : keyList) {
            Double value = pairToValue.get(key);
            System.out.println(key + "\t" + value);
        }
        System.out.println("\nComponents between GSE4922 to nejm:");
        pairToValue = calculatePValuesForTwoComponentSets(gse4922Comps, 
                                                          nejmComps);
        keyList = new ArrayList<String>(pairToValue.keySet());
        Collections.sort(keyList);
        keyList = new ArrayList<String>(pairToValue.keySet());
        Collections.sort(keyList);
        for (String key : keyList) {
            Double value = pairToValue.get(key);
            System.out.println(key + "\t" + value);
        }
    }
    
    @Test
    public void checkWithKnownGenes() throws Exception {
        Set<String> knownGenes = getChuangKnownBreastGenes();
        List<GraphComponent> gse2034Components = loadGSE2034Comonents();
        List<GraphComponent> gse4922Components = loadGSE4922Components();
        List<GraphComponent> nejmComponents = loadNejmComponents();
        System.out.println("Check how many components have been picked up:");
        System.out.println("Total compnents in GSE2034: " + gse2034Components.size());
        System.out.println("Total components in GSE4922: " + gse4922Components.size());
        System.out.println("Total components in van de Vijver: " + nejmComponents.size());
        // Get all ids in the components
        Set<String> idsInGse2034 = getIdsFromComponents(gse2034Components);
        Set<String> idsInGse4922 = getIdsFromComponents(gse4922Components);
        Set<String> idsInNejm = getIdsFromComponents(nejmComponents);
        // Check the total ids in components:
        System.out.println("Total ids in all components:");
        System.out.println("Ids in GSE2034: " + idsInGse2034.size());
        System.out.println("Ids in GSE4922: " + idsInGse4922.size());
        System.out.println("Ids in van de Vijver: " + idsInNejm.size());
        // Check shared ids in two data sets
        Set<String> copy2034 = new HashSet<String>(idsInGse2034);
        Set<String> copy4922 = new HashSet<String>(idsInGse4922);
        Set<String> copyNejm = new HashSet<String>(idsInNejm);
        // Check shared genes
        copy2034.retainAll(knownGenes);
        copy4922.retainAll(knownGenes);
        copyNejm.retainAll(knownGenes);
        System.out.println("Check coverage of known genes:");
        System.out.println("In GSE2034 components: " + copy2034.size());
        System.out.println("In GSE4922 components: " + copy4922.size());
        System.out.println("In van de Vijver: " + copyNejm.size());
        copy2034 = new HashSet<String>(idsInGse2034);
        copy2034.retainAll(idsInGse4922);
        System.out.println("Shared in GSE2034 and GSE4922: " + copy2034.size());
        copy2034.retainAll(knownGenes);
        System.out.println("In GSE2034 and GSE4922 components: " + copy2034.size());
        copy2034 = new HashSet<String>(idsInGse2034);
        copy2034.retainAll(idsInGse4922);
        copy2034.retainAll(idsInNejm);
        int allShared = copy2034.size();
        System.out.println("Shared in all three: " + copy2034.size());
        copy2034.retainAll(knownGenes);
        int hit = copy2034.size();
        System.out.println("In All three: " + copy2034.size());
        int totalNetworkSize = 9815;
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalNetworkSize,
                                                                    knownGenes.size(),
                                                                    allShared,
                                                                    hit);
        System.out.println("p value: " + pvalue);
    }

    @Test
    public void filterGSE4922BasedOnSampleInfo() throws IOException {
        Map<String, String> sampleToType = dataAnalyzer.loadGSE4922PatientInfo();
        System.out.println("SampleToType: " + sampleToType.size());
        String inFileName = dirName + "GSE4922Filtered.txt";
        String outFileName = dirName + "GSE4922FilteredOnSamples.txt";
        fu.setInput(inFileName);
        String line = fu.readLine();
        List<String> samples = new ArrayList<String>();
        dataAnalyzer.extractGSE2034Samples(line, samples);
        System.out.println("Total sample: " + samples.size());
        int lastIndex = 0;
        for (int i = 0; i < samples.size(); i++) {
            String sample = samples.get(i);
            if (sampleToType.containsKey(sample))
                lastIndex = i;
        }
        System.out.println("Last Index: " + lastIndex);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("gene");
        for (int i = 0; i <= lastIndex; i++)
            builder.append("\t").append(samples.get(i));
        outFu.printLine(builder.toString());
        lastIndex ++; // Need to consider the first gene element
        while ((line = fu.readLine()) != null) {
            builder.setLength(0);
            String[] tokens = line.split("\t");
            builder.append(tokens[0]);
            List<String> values = new ArrayList<String>();
            for (int i = 1; i <= lastIndex; i++) {
                builder.append("\t").append(tokens[i]);
            }
            outFu.printLine(builder.toString());
        }
        fu.close();
        outFu.close();
    }
    
    private Map<String, Double> calculatePValuesForTwoComponentSets(List<GraphComponent> comps1,
                                                                    List<GraphComponent> comps2) throws MathException {
        Map<String, Double> pairToPValue = new HashMap<String, Double>();
        for (int i = 0; i < comps1.size(); i++) {
            GraphComponent comp1 = comps1.get(i);
            Set<String> ids1 = comp1.getAllNodes();
            for (int j = 0; j < comps2.size(); j++) {
                GraphComponent comp2 = comps2.get(j);
                Set<String> ids2 = comp2.getAllNodes();
                int shared = sharedElements(ids1, ids2);
                double p = MathUtilities.calculateHypergeometricPValue(9500, 
                                                                       ids1.size(), 
                                                                       ids2.size(), 
                                                                       shared);
                if (p > MAX_P_VALUE)
                    continue;
                String key = i + " " + j;
                pairToPValue.put(key, p);
            }
        }
        return pairToPValue;
    }
    
    /**
     * This method is used to analyze overlapping among two clusters.
     * @throws Exception
     */
    @Test
    public void analyzeClusterOverlappingWithPValue() throws Exception {
        // Get the total genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        int total = fiGenes.size();
        String fileName = dirName + "FIsWithWeightFromTTestMCLClustersGSE2034_I80_091509.txt";
        MCLResultsAnalyzer mclAnalyzer = new MCLResultsAnalyzer();
        List<Set<String>> gse2034Clusters = mclAnalyzer.loadMCLClusters(fileName);
        Map<String, Double> fiToWeight1 = loadFIsWithWeights(dirName + "FIsWithWeightFromTTestGSE2034_091509.txt");
        fileName = dirName + "FIsWithWeightFromTTestMCLClustersGSE4922_I80_091509.txt";
        // Generate with random test
        //fileName = dirName + "FIsWithRandomWeightsMCLClustersGSE4922_I80_091409.txt";
        List<Set<String>> gse4922Clusters = mclAnalyzer.loadMCLClusters(fileName);
        Map<String, Double> fiToWeight2 = loadFIsWithWeights(dirName + "FIsWithWeightFromTTestGSE4922_091509.txt");
        System.out.println("GSE2034 Cluster\tSize\tAverage_Weight\tGSE 4922 Cluster\tSize\tAverage_Weight\tShared\tP-value");
        // Want to calculate pair-wise cluster comparison
        for (int i = 0; i < gse2034Clusters.size(); i++) {
            Set<String> cluster1 = gse2034Clusters.get(i);
            if (cluster1.size() < 5)
                continue;
            double averageWeight1 = calculateAverageWeight(fiToWeight1, cluster1, fis);
            if (averageWeight1 < 0.30)
                continue; // Based on a random test
            for (int j = 0; j < gse4922Clusters.size(); j++) {
                Set<String> cluster2 = gse4922Clusters.get(j);
                if (cluster2.size() < 5)
                    continue;
                double averageWeight2 = calculateAverageWeight(fiToWeight2, cluster2, fis);
                // Don't apply the score filter here
                //if (averageWeight2 < 0.30)
                //    continue;
                Set<String> shared = new HashSet<String>(cluster2);
                shared.retainAll(cluster1);
                double pvalue = MathUtilities.calculateHypergeometricPValue(total, 
                                                                            cluster1.size(),
                                                                            cluster2.size(), 
                                                                            shared.size());
                if (pvalue > 1.0E-6) // Use this cutoff for Benferonni test
                    continue;
                System.out.println(i + "\t" + 
                                   cluster1.size() + "\t" +
                                   averageWeight1 + "\t" + 
                                   j + "\t" + 
                                   cluster2.size() + "\t" +
                                   averageWeight2 + "\t" +
                                   shared.size() + "\t" +
                                   pvalue);
            }
        }
    }
    
    private double calculateAverageWeight(Map<String, Double> fiToWeight,
                                          Set<String> cluster,
                                          Set<String> fis) {
        Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
        if (fisInCluster.size() == 0)
            return 0.0;
        double total = 0.0d;
        for (String fi : fisInCluster) {
            double weight = fiToWeight.get(fi);
            total += weight;
        }
        return total / fisInCluster.size();
    }
    
    private Map<String, Double> loadFIsWithWeights(String fileName) throws IOException {
        fu.setInput(fileName);
        Map<String, Double> fiToWeight = new HashMap<String, Double>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fiToWeight.put(tokens[0] + "\t" + tokens[1],
                           new Double(tokens[2]));
        }
        fu.close();
        return fiToWeight;
    }
    
    @Test
    public void analyzeComponentsOverlapping() throws IOException {
        // Want to check the total size of the graph
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        int maxSize = (int) (totalGenes.size() * 0.05); // Pick up at most 0.05
        // From Wang et al
        String compFileName = dirName + "Components2_1_Norm.txt";
        // Get the value from permutation test
        double minimumMI = 0.093257329978305d;
        List<GraphComponent> gse2034Components = pickUpComponents(compFileName, 
                                                                  minimumMI, 
                                                                  maxSize);
        // For GSE4922
        compFileName = dirName + "Components2_1_GSE4922_Norm.txt";
        minimumMI = 0.09405237106086342;
        List<GraphComponent> gse4922Components = pickUpComponents(compFileName,
                                                                  minimumMI,
                                                                  maxSize);
        // For van de Vijver data set
        compFileName = dirName + "Components2_1_Nejm.txt";
        minimumMI = 0.11537420328179235d;
        List<GraphComponent> nejmComponents = pickUpComponents(compFileName,
                                                               minimumMI,
                                                               maxSize);
        System.out.println("Check how many components have been picked up:");
        System.out.println("Total compnents in GSE2034: " + gse2034Components.size());
        System.out.println("Total components in GSE4922: " + gse4922Components.size());
        System.out.println("Total components in van de Vijver: " + nejmComponents.size());
        // Merging
        List<GraphComponent> gse2034Merged = mergeComponents(gse2034Components);
        List<GraphComponent> gse4922Merged = mergeComponents(gse4922Components);
        List<GraphComponent> nejmMerged = mergeComponents(nejmComponents);
        System.out.println("After component merging:");
        System.out.println("Total components in GSE2034: " + gse2034Merged.size());
        System.out.println("Total components in GSE4922: " + gse4922Merged.size());
        System.out.println("Total components in van de Vijver: " + nejmMerged.size());
        System.out.println("Check component size:");
        System.out.println("Components in GSE2034:");
        checkComponentSize(gse2034Merged);
        System.out.println("Components in GSE4922:");
        checkComponentSize(gse4922Merged);
        System.out.println("Components in van de Vijver:");
        checkComponentSize(nejmMerged);
        // Get all ids in the components
        Set<String> idsInGse2034 = getIdsFromComponents(gse2034Merged);
        Set<String> idsInGse4922 = getIdsFromComponents(gse4922Merged);
        Set<String> idsInNejm = getIdsFromComponents(nejmMerged);
        // Check the total ids in components:
        System.out.println("Total ids in all components:");
        System.out.println("Ids in GSE2034: " + idsInGse2034.size());
        System.out.println("Ids in GSE4022: " + idsInGse4922.size());
        System.out.println("Ids in van de Vijver: " + idsInNejm.size());
        // Check shared ids in two data sets
        Set<String> copy2034 = new HashSet<String>(idsInGse2034);
        Set<String> copy4922 = new HashSet<String>(idsInGse4922);
        Set<String> copyNejm = new HashSet<String>(idsInNejm);
        copy2034.retainAll(copy4922);
        System.out.println("Shared genes between GSE2034 and GSE4922: " + copy2034.size());
        copy2034 = new HashSet<String>(idsInGse2034);
        copy2034.retainAll(idsInNejm);
        System.out.println("Shared genes between GSE2034 and van de Vijver: " + copy2034.size());
        copy4922.retainAll(idsInNejm);
        System.out.println("Shared genes between GSE4922 and van de Vijver: " + copy4922.size());
        // All shared
        copy4922.retainAll(idsInGse2034);
        System.out.println("Genes shared in all three data sets: " + copy4922.size());
        System.out.println("Shared genes are: " + copy4922);
        for (GraphComponent comp1 : gse2034Merged) {
            for (int i = 0; i < nejmMerged.size(); i++) {
                GraphComponent nejmComp = nejmMerged.get(i);
                if (interact(comp1, nejmComp, fis)) {
                    System.out.println("Nejm component " + i + " can interact with gse2034.");
                }
            }
        }
        for (GraphComponent comp1 : gse4922Merged) {
            for (int i = 0; i < nejmMerged.size(); i++) {
                GraphComponent nejmComp = nejmMerged.get(i);
                if (interact(comp1, nejmComp, fis)) {
                    System.out.println("Nejm component " + i + " can interact with gse4922.");
                }
            }
        }
    }
    
    private void checkComponentSize(List<GraphComponent> comps) {
        for (int i = 0; i < comps.size(); i++) {
            GraphComponent comp = comps.get(i);
            System.out.println("Comp " +  i + ": " + comp.getAllNodes().size());
        }
    }
    
    private Set<String> getIdsFromComponents(List<GraphComponent> comps) {
        Set<String> ids = new HashSet<String>();
        for (GraphComponent comp : comps) {
            Set<String> compIds = comp.getAllNodes();
            ids.addAll(compIds);
        }
        return ids;
    }
    
    private List<GraphComponent> mergeComponents(List<GraphComponent> comps) throws IOException {
        int preSize = comps.size();
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        while (true) {
            GraphComponent target1 = null, target2 = null;
            for (int i = 0; i < comps.size() - 1; i++) {
                GraphComponent comp1 = comps.get(i);
                target2 = null;
                for (int j = i + 1; j < comps.size(); j++) {
                    GraphComponent comp2 = comps.get(j);
                    if (interact(comp1, comp2, fis)) {
                        target2 = comp2;
                        break;
                    }
                }
                if (target2 != null) {
                    target1 = comp1;
                    break;
                }
            }
            if (target2 == null)
                break;
            // Merging
            Set<String> nodes2 = target2.getAllNodes();
            for (String node : nodes2)
                target1.addNode(node);
            comps.remove(target2);
        }
        return comps;
    }
    
    /**
     * Used to annotate gene sets.
     * @throws Exception
     */
    @Test
    public void annotateGeneSet() throws Exception {
//        Set<String> genes = new HashSet<String>();
//        String fileName= dirName + "GeneSharedInThreeDatasets.txt";
//        fu.setInput(fileName);
//        String line = fu.readLine();
//        fu.close();
//        String[] tokens = line.split(",");
//        for (String token : tokens)
//            genes.add(token.trim());
        //Set<String> genes = getKnownBreastCancerGenes();
        // From Wang et al
        List<GraphComponent> gse2034Components = loadGSE2034Comonents();
        // For GSE4922
        List<GraphComponent> gse4922Components = loadGSE4922Components();
        // For van de Vijver data set
        List<GraphComponent> nejmComponents = loadNejmComponents();
        int index = 0;
        Set<String> allPathways = new HashSet<String>();
        for (GraphComponent comp : nejmComponents) {
            if (index != 50 && index != 55) {
                index ++;
                continue;
            }
            System.out.println("Component " + index);
            index ++;
            Set<String> genes = comp.getAllNodes();
            // Need to map to protein ids
            UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
            Map<String, String> geneNameToProteinIds = uniAnalyzer.mapGeneNamesToUniProt(genes);
            System.out.println("Total genes: " + genes.size());
            System.out.println("Total proteins: " + geneNameToProteinIds.values().size());
            PathwayBasedAnnotator topicAnalyzer = new PathwayBasedAnnotator();
//            topicAnalyzer.setHopNumber(0);
            topicAnalyzer.setPValueThreshold(MAX_P_VALUE);
            //FileOutputStream fos = new FileOutputStream(dirName + "tmp.txt");
            //PrintStream ps = new PrintStream(fos);
            PrintStream ps = System.out;
            topicAnalyzer.annotateGenesWithFDR(geneNameToProteinIds.values(), AnnotationType.Pathway);
            //ps.close();
            //fos.close();
            //Set<String> pathways =  parsePathways(dirName + "tmp.txt");
            //allPathways.addAll(pathways);
        }
        System.out.println("Pathways: " + allPathways.size());
        for (String pathway : allPathways)
            System.out.println(pathway);
    }

    private void loadTotalGenes() throws IOException {
        // Want to check the total size of the graph
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        allFIs = fu.loadInteractions(intFileName);
        totalGenes = InteractionUtilities.grepIDsFromInteractions(allFIs);
    }
    
    private int sharedElements(Set<String> set1,
                               Set<String> set2) {
        int shared = 0;
        for (String element1 : set1) {
            if (set2.contains(element1))
                shared ++;
        }
        return shared;
    }
    
    /**
     * Check component overlapping with known breast cancer genes.
     * @throws Exception
     */
    @Test
    public void checkCompFIsWithKnownGenes() throws Exception {
        if (totalGenes == null)
            loadTotalGenes();
        Set<String> chuangKnownGenes = getChuangKnownBreastGenes();
        Set<String> sjobormGenes = getSjoblomKnownBreastCancerGenes();
        // shared genes
        int shared = sharedElements(chuangKnownGenes, sjobormGenes);
        System.out.println("Shared known genes: " + shared);
        // Want to check the total size of the graph
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        List<GraphComponent> gse2034Comps = loadGSE2034Comonents();
        List<GraphComponent> gse4922Comps = loadGSE4922Components();
        List<GraphComponent> nejmComps = loadNejmComponents();
        System.out.println("Checked with chuang genes:");
        checkFIsWithKnownGenes(chuangKnownGenes,
                               totalGenes,
                               gse4922Comps);
        checkFIsWithKnownGenes(chuangKnownGenes,
                               totalGenes, 
                               gse4922Comps, 
                               fis,
                               true);
        System.out.println("\nChecked with sjoblom genes:");
        checkFIsWithKnownGenes(sjobormGenes, 
                               totalGenes, 
                               gse4922Comps);
        checkFIsWithKnownGenes(sjobormGenes, 
                               totalGenes, 
                               gse4922Comps, 
                               fis,
                               true);
    }
    
    private void checkFIsWithKnownGenes(Set<String> chuangKnownGenes,
                                        Set<String> totalGenes,
                                        List<GraphComponent> comps)
            throws MathException {
        System.out.println("Simple gene overlapping");
        // Check overlapping with known breast cancer genes
        System.out.println("Component\tTotalGenesInComponent\tCompGenesInKnownGenes\tPvalue\tHitGenes");
        for (int i = 0; i < comps.size(); i++) {
            GraphComponent comp = comps.get(i);
            Set<String> nodes = comp.getAllNodes();
            Set<String> copy = new HashSet<String>(nodes);
            copy.retainAll(chuangKnownGenes);
            Double pValue = null;
            if (copy.size() > 0) {
                pValue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(), 
                                                                     nodes.size(), 
                                                                     chuangKnownGenes.size(), 
                                                                     copy.size());
            }
            else
                pValue = 1.0d;
            System.out.println("Component " + i + "\t" + nodes.size() + "\t" +
                               copy.size() + 
                               (pValue == null ? "" : ("\t" + pValue)) + "\t" +
                               copy);
        }
    }
    
    /**
     * Random test a set of genes in GraphComponents and their FI partners.
     * @throws Exception
     */
    @Test
    public void randomCheckCompFIsWithKnownGenes() throws Exception {
        if (totalGenes == null)
            loadTotalGenes();
        Set<String> chuangKnownGenes = getChuangKnownBreastGenes();
        // Want to check the total size of the graph
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
//        List<GraphComponent> gse2034Comps = loadGSE2034Comonents();
//        String pvalueFileName = dirName + "RandomTest100FICompPValuesGSE2034.txt";
//        List<GraphComponent> comps = gse2034Comps;
//        List<GraphComponent> gse4922Comps = loadGSE4922Components();
//        String pvalueFileName = dirName + "RandomTest100FICompPValuesGSE4922.txt";
//        List<GraphComponent> comps = gse4922Comps;
        List<GraphComponent> nejmComps = loadNejmComponents();
        String pvalueFileName = dirName + "RandomTest1000FICompPValuesNejm.txt";
        List<GraphComponent> comps = nejmComps;
        long time1 = System.currentTimeMillis();
        List<Double> allPValues = new ArrayList<Double>();
        int permutation = 1000;
        for (int i = 0; i < permutation; i++) {
            // Random pick up genes
            Set<String> randomGenes = sampleList(new ArrayList<String>(totalGenes),
                                                 chuangKnownGenes.size());
//          checkFIsWithKnownGenes(randomGenes,
//          totalGenes,
//          gse2034Comps);
            List<Double> pvalues = checkFIsWithKnownGenes(randomGenes,
                                                          totalGenes, 
                                                          comps, 
                                                          fis,
                                                          false);
            // Print out the smallest p-value
            System.out.println(i + " minum: " + pvalues.get(0));
            allPValues.addAll(pvalues);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for one run: " + (time2 - time1));
        Collections.sort(allPValues);
        fu.setOutput(pvalueFileName);
        for (Double pValue : allPValues)
            fu.printLine(pValue.toString());
        fu.close();
    }
    
    /**
     * Calculate p values for a list of GraphComponents against FIs.
     * @param knownGenes
     * @param totalGenes
     * @param comps
     * @param fis
     * @return
     * @throws Exception
     */
    private List<Double> checkFIsWithKnownGenes(Set<String> knownGenes,
                                        Set<String> totalGenes,
                                        List<GraphComponent> comps,
                                        Set<String> fis,
                                        boolean needOutput) throws Exception {
        if (needOutput) {
            System.out.println("based on FIs:");
            System.out.println("Component\tTotalFIPatners\tPartnerInKnownGenes\tPvalue\tHitGenes");
        }
        // Check how many knowns genes interact with comps
        List<Double> pValues = new ArrayList<Double>();
        for (int i = 0; i < comps.size(); i++) {
            GraphComponent comp = (GraphComponent) comps.get(i);
            Set<String> nodes = comp.getAllNodes();
            Set<String> sharedOrFiKnown = new HashSet<String>();
            for (String known : knownGenes) {
                if (!nodes.contains(known) && interact(known, nodes, fis))
                    sharedOrFiKnown.add(known);
            }
            Set<String> totalFiPartners = getFIInteractors(nodes, fis);
            double pValue = 1.0d;
            if (sharedOrFiKnown.size() > 0)
                pValue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(), 
                                                                     totalFiPartners.size(), 
                                                                     knownGenes.size(),
                                                                     sharedOrFiKnown.size());
            if (needOutput) {
                System.out.println("Component " + i + "\t" + totalFiPartners.size() + "\t" +
                                   sharedOrFiKnown.size() +  "\t" +
                                   pValue + "\t" +
                                   sharedOrFiKnown);
            }
            pValues.add(pValue);
        }
        Collections.sort(pValues);
        return pValues;
    }
    
    
    private Set<String> getFIInteractors(Set<String> ids,
                                         Set<String> fis) {
        Set<String> interactors = new HashSet<String>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            if (ids.contains(id1) || ids.contains(id2)) {
                interactors.add(id1);
                interactors.add(id2);
            }
        }
        interactors.removeAll(ids);
        return interactors;
    }
    
    private Set<String> getFIInteractors(String id,
                                         Set<String> fis) {
        Set<String> interactors = new HashSet<String>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            if (id1.equals(id))
                interactors.add(id2);
            else if (id2.equals(id))
                interactors.add(id1);
        }
        return interactors;
    }
    
    private Set<String> getChuangKnownBreastGenes() throws IOException {
        // Load known genes
        String fileName = dirName + "ListOfBreastCancerGenesFromChuang.txt";
        String delimit = "\t";
        return getKnownBreastCancerGenes(fileName, delimit);
    }
    
    private Set<String> getSjoblomKnownBreastCancerGenes() throws IOException {
        String fileName = dirName + "BreastCancerGeneFromSjoblom.txt";
        String delimit = " ";
        return getKnownBreastCancerGenes(fileName, delimit);
    }

    /**
     * The returned genes have been filtered based on the biggest component used for
     * component searching.
     * @param fileName
     * @param delimit
     * @return
     * @throws IOException
     */
    private Set<String> getKnownBreastCancerGenes(String fileName,
                                                  String delimit)
            throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> knownGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split(delimit);
            String gene = tokens[0];
            if (gene.contains("(")) {
                int index = gene.indexOf("(");
                knownGenes.add(gene.substring(0, index));
                knownGenes.add(gene.substring(0, gene.length() - 1));
            }
            else
                knownGenes.add(gene);
        }
        fu.close();
        System.out.println("Total known genes: " + knownGenes.size());
        EntrezGeneAnalyzer geneAnalyzer = new EntrezGeneAnalyzer();
        knownGenes = geneAnalyzer.normalizeGeneNames(knownGenes);
        System.out.println("After name normalization: " + knownGenes.size());
        // Need to do another filtering
        Set<String> fis = fu.loadInteractions(R3Constants.RESULT_DIR + "FI73InGene_061008.txt");
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total size of network: " + allGenes.size());
        knownGenes.retainAll(allGenes);
        System.out.println("After network filtering: " + knownGenes.size());
        return knownGenes;
    }
    
    private List<GraphComponent> loadGSE2034Comonents() throws IOException {
        if (totalGenes == null)
            loadTotalGenes();
        int maxSize = (int)(totalGenes.size() * 0.05);
        // From Wang et al
        String compFileName = dirName + "Components1_1_GSE2023.txt";
        //String compFileName = dirName + "Components2_1_Norm.txt";
        //String compFileName = dirName + "Components2_1_GSE2023_ID_Permutated_1.txt";
        // Get the value from perumation test
        double minimumMI = 0.093257329978305d;
        List<GraphComponent> gse2034Components = pickUpComponents(compFileName, 
                                                                  minimumMI, 
                                                                  maxSize);
        return gse2034Components;
    }
    
    private List<GraphComponent> loadGSE4922Components() throws IOException {
        if (totalGenes == null)
            loadTotalGenes();
        int maxSize = (int)(totalGenes.size() * 0.05);
        // For GSE4922
        String compFileName = dirName + "Components2_1_GSE4922_Norm.txt";
        //String compFileName = dirName + "Components2_1_GSE4922_ID_Permutated.txt";
        double minimumMI = 0.09405237106086342;
        List<GraphComponent> gse4922Components = pickUpComponents(compFileName,
                                                                  minimumMI,
                                                                  maxSize);
        return gse4922Components;
    }
    
    private List<GraphComponent> loadNejmComponents() throws IOException {
        if (totalGenes == null)
            loadTotalGenes();
        int maxSize = (int)(totalGenes.size() * 0.05);
        // For van de Vijver data set
        String compFileName = dirName + "Components2_1_Nejm.txt";
        double minimumMI = 0.11537420328179235d;
        List<GraphComponent> nejmComponents = pickUpComponents(compFileName,
                                                               minimumMI,
                                                               maxSize);
        return nejmComponents;
    }
    
    
    private boolean interact(GraphComponent comp1,
                             GraphComponent comp2,
                             Set<String> fis) {
        // If shared, two components must interact
        if (share(comp1, comp2, 0.0d))
            return true;
        int compare = 0;
        String fi = null;
        for (String id1 : comp1.getAllNodes()) {
            for (String id2 : comp2.getAllNodes()) {
                compare = id1.compareTo(id2);
                if (compare < 0) {
                    fi = id1 + "\t" + id2;
                }
                else if (compare > 0)
                    fi = id2 + "\t" + id1;
                if (fi != null &&
                    fis.contains(fi))
                    return true;
            }
        }
        return false;
    }
    
    
    private boolean share(GraphComponent comp1,
                          GraphComponent comp2,
                          double cutoff) {
        Set<String> nodes1 = comp1.getAllNodes();
        Set<String> nodes2 = comp2.getAllNodes();
        Set<String> copy1 = new HashSet<String>(nodes1);
        Set<String> copy2 = new HashSet<String>(nodes2);
        copy1.retainAll(copy2);
        int smaller = Math.min(nodes1.size(), nodes2.size());
        double shared = (double) copy1.size() / smaller;
        return shared > cutoff;
    }
    
    private boolean interact(String target,
                             Set<String> ids,
                             Set<String> fis) {
        for (String id : ids) {
            int compare = target.compareTo(id);
            String fi = null;
            if (compare < 0)
                fi = target + "\t" + id;
            else
                fi = id + "\t" + target;
            if (fis.contains(fi))
                return true;
        }
        return false;
    }
    
    
    /**
     * This helper method is used to pick up GraphComponents so that the 
     * total GraphComponents will have the required idSize (should not exceeding),
     * and their MI values should higher than minimumMI. The second condition is
     * most important.
     * @param components
     * @param minimumMi
     * @param idSize
     * @return
     */
    public List<GraphComponent> pickUpComponents(String compFileName,
                                                 double minimumMi,
                                                 int idSize) throws IOException {
        List<GraphComponent> components = loadComponentsInOrder(compFileName);
        List<GraphComponent> filtered = new ArrayList<GraphComponent>();
        Set<String> totalIds = new HashSet<String>();
        for (int i = 0; i < components.size(); i++) {
            GraphComponent comp = components.get(i);
            // Check the score
            if (comp.getScore() < minimumMi)
                break;
            // Check total ids
            // Do a peek
            Set<String> copy = new HashSet<String>(totalIds);
            copy.addAll(comp.getAllNodes());
            if (copy.size() > idSize)
                break;
            totalIds.addAll(comp.getAllNodes());
            filtered.add(comp);
        }
        return filtered;
    }
    
    /**
     * Load a list of GraphComponents sorted by scores.
     * @return
     * @throws IOException
     */
    public List<GraphComponent> loadComponentsInOrder(String fileName) throws IOException {
        List<GraphComponent> components = new ArrayList<GraphComponent>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            GraphComponent component = new GraphComponent();
            component.setId(new Integer(tokens[0]));
            for (int i = 3; i < tokens.length; i++)
                component.addNode(tokens[i]);
            component.setScore(new Double(tokens[1]));
            components.add(component);
        }
        fu.close();
        // Need to sort 
        Collections.sort(components, new Comparator<GraphComponent>() {
            public int compare(GraphComponent comp1,
                               GraphComponent comp2) {
                Double score1 = comp1.getScore();
                Double score2 = comp2.getScore();
                return score2.compareTo(score1);
            }
        });
        return components;
    }

    private List<GraphComponent> filterComponents(double miCutOff,
                                                  List<GraphComponent> components) {
        List<GraphComponent> filteredComps = new ArrayList<GraphComponent>();
        for (GraphComponent comp : components) {
            if (comp.getScore() >= miCutOff &&
                comp.getAllNodes().size() >= CUT_OFF_OF_COMP_SIZE)
                filteredComps.add(comp);
        }
        return filteredComps;
    }

    public void filterRedundantComponents(List<GraphComponent> components) {
        // If a lower score component has genes contained by higher score component,
        // don't use this component
        List<GraphComponent> toBeRemoved = new ArrayList<GraphComponent>();
        for (int i = components.size() - 1; i > 0; i--) {
            GraphComponent comp = components.get(i);
            for (int j = i - 1; j >= 0; j--) {
                GraphComponent better = components.get(j);
                if (share(better, comp, 0.75)) {
                    toBeRemoved.add(comp);
                }
            }
        }
        components.removeAll(toBeRemoved);
    }
    
    /**
     * This method is used to permutate test samples in a components for p-values calculation.
     * @throws Exception
     */
    @Test
    public void permutateSamplesInComponents() throws Exception {
        long time1 = System.currentTimeMillis();
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        String outFileName = dirName + "Components2_1_GSE2023_NO_P_I_MI_random_sample_1000.txt";
        String compFileName = dirName + "Components2_1_GSE2023_NO_P_I.txt";
//        // For van de Vijver data set
//        miCalculator.setSampleToPhenotypeInfo(loadNejmPatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
//        String outFileName = dirName + "Components2_1_Nejm_MI_random_sample_1000.txt";
//        String compFileName = dirName + "Components2_1_Nejm.txt";
        // GSE4922
//        miCalculator.setSampleToPhenotypeInfo(loadGSE4922PatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore.txt");
//        String outFileName = dirName + "Components2_1_GSE4922_Norm_MI_random_sample_1000.txt";
//        String compFileName = dirName + "Components2_1_GSE4922_Norm.txt";
        
        Map<Integer, Set<String>> clusterToGenes = getComponentToGenes(compFileName);
        Map<Integer, List<String>> clusterToDiscValues = new HashMap<Integer, List<String>>();
        for (Iterator<Integer> it = clusterToGenes.keySet().iterator(); it.hasNext();) {
            Integer cluster = it.next();
            Set<String> genes = clusterToGenes.get(cluster);
            if (genes.size() < CUT_OFF_OF_COMP_SIZE)
                continue;
            GraphComponent component = new GraphComponent();
            for (String gene : genes)
                component.addNode(gene);
            List<String> disValues = miCalculator.calculateDiscValues(component);
            clusterToDiscValues.put(cluster, disValues);
        }
        int numberOfPermutation = 1000;
        // Try to use a randomized sample list
        List<String> originalSamples = miCalculator.getSamples();
        List<Double> scores = new ArrayList<Double>();
        for (int i = 0; i < numberOfPermutation; i++) {
            long time3 = System.currentTimeMillis();
            List<String> samples = permutateList(originalSamples);
            miCalculator.setSamples(samples);
            for (Iterator<Integer> it = clusterToDiscValues.keySet().iterator(); 
                 it.hasNext();) {
                Integer cluster = it.next();
                List<String> discValues = clusterToDiscValues.get(cluster);
                double mi = miCalculator.calculateScore(discValues);
                scores.add(mi);
            }
            long time4 = System.currentTimeMillis();
            System.out.println("Time for loop " + i + ": " + (time4 - time3));
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
        fu.setOutput(outFileName);
        Collections.sort(scores, new Comparator<Double>() {
            public int compare(Double score1, Double score2) {
                return score2.compareTo(score1);
            }
        });
        for (Double score : scores)
            fu.printLine(score.toString());
        fu.close();
    }
    
    private Map<Integer, Set<String>> getComponentToGenes(String fileName) throws IOException {
        Map<Integer, Set<String>> compToGenes = new HashMap<Integer, Set<String>>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Set<String> genes = new HashSet<String>();
            for (int i = 2; i < tokens.length; i++)
                genes.add(tokens[i]);
            compToGenes.put(Integer.parseInt(tokens[0]),
                            genes);
        }
        fu.close();
        return compToGenes;
    }
    
    /**
     * Permutate a list (i.e. reshuffle the elements in the list).
     * @param samples
     * @return
     */
    private List<String> permutateList(List<String> samples) {
        ObjectArrayList list = new ObjectArrayList();
        list.addAllOf(samples);
        list.shuffle();
        List<String> rtn = new ArrayList<String>(samples.size());
        for (int i = 0; i < list.size(); i++)
            rtn.add(list.get(i).toString());
        return rtn;
    }
    
    /**
     * Sample from a list for a pre-defined size.
     * @param list
     * @param sampleSize
     * @return
     */
    private Set<String> sampleList(List<String> list,
                                   int sampleSize) {
        int size = list.size();
        Set<String> sample = new HashSet<String>();
        int index = 0;
        while (sample.size() < sampleSize) {
            index = (int) (Math.random() * size);
            sample.add(list.get(index));
        }
        return sample;
    }
    
    @Test
    public void randomCheckComponentMerging() throws IOException {
        // Want to check the total size of the graph
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        int maxSize = (int) (totalGenes.size() * 0.01); // Pick up at most 0.05

        // From Wang et al
        String compFileName = dirName + "Components2_1_Norm.txt";
        Set<GraphComponent> gse2034Components = randomPickupComponents(compFileName, 
                                                                        maxSize);
        // For GSE4922
        compFileName = dirName + "Components2_1_GSE4922_Norm.txt";
        Set<GraphComponent> gse4922Components = randomPickupComponents(compFileName,
                                                                  maxSize);
        // For van de Vijver data set
        compFileName = dirName + "Components2_1_Nejm.txt";
        Set<GraphComponent> nejmComponents = randomPickupComponents(compFileName,
                                                               maxSize);
        System.out.println("Check how many components have been picked up:");
        System.out.println("Total compnents in GSE2034: " + gse2034Components.size());
        System.out.println("Total components in GSE4922: " + gse4922Components.size());
        System.out.println("Total components in van de Vijver: " + nejmComponents.size());
        // Merging
        List<GraphComponent> gse2034Merged = mergeComponents(new ArrayList<GraphComponent>(gse2034Components));
        List<GraphComponent> gse4922Merged = mergeComponents(new ArrayList<GraphComponent>(gse4922Components));
        List<GraphComponent> nejmMerged = mergeComponents(new ArrayList<GraphComponent>(nejmComponents));
        System.out.println("After component merging:");
        System.out.println("Total components in GSE2034: " + gse2034Merged.size());
        System.out.println("Total components in GSE4922: " + gse4922Merged.size());
        System.out.println("Total components in van de Vijver: " + nejmMerged.size());
        System.out.println("Check component size:");
        System.out.println("Components in GSE2034:");
        checkComponentSize(gse2034Merged);
        System.out.println("Components in GSE4922:");
        checkComponentSize(gse4922Merged);
        System.out.println("Components in van de Vijver:");
        checkComponentSize(nejmMerged);
    }
    
    @Test
    public void anlyzePathwayOverlapping() throws IOException {
        String fileName = dirName + "HitPathwayComparision.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> gse2043Pathways = new HashSet<String>();
        Set<String> gse4922Pathways = new HashSet<String>();
        Set<String> njemPathways = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].length() > 0)
                gse2043Pathways.add(tokens[0]);
            if (tokens.length > 1 && tokens[1].length() > 0)
                gse4922Pathways.add(tokens[1]);
            if (tokens.length > 2 && tokens[2].length() > 0)
                njemPathways.add(tokens[2]);
        }
        System.out.println("Pathways in gse2034: " + gse2043Pathways.size());
        System.out.println("Pathways in gse4922: " + gse4922Pathways.size());
        System.out.println("Pathways in nejm: " + njemPathways.size());
        Set<String> gse2034Copy = new HashSet<String>(gse2043Pathways);
        gse2034Copy.retainAll(gse4922Pathways);
        System.out.println("shared in gse2034 and 4922: " + gse2034Copy.size());
        gse2034Copy = new HashSet<String>(gse2043Pathways);
        gse2034Copy.retainAll(njemPathways);
        System.out.println("shared in gse2034 and nejm: " + gse2034Copy.size());
        Set<String> gse4922Copy = new HashSet<String>(gse4922Pathways);
        gse4922Copy.retainAll(njemPathways);
        System.out.println("shared in gse4922 and nejm: " + gse4922Copy.size());
        gse4922Copy.retainAll(gse2043Pathways);
        System.out.println("shared in all: " + gse4922Copy.size());
    }
    
    private Set<String> parsePathways(String tmpFileName) throws IOException {
        Set<String> pathways = new HashSet<String>();
        fu.setInput(tmpFileName);
        String line = null;
        boolean isInData = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Pathway")) {
                isInData = true;
            }
            else if (line.startsWith("Total ids:")) {
                isInData = false;
            }
            else if (isInData) {
                String[] tokens = line.split("\t");
                pathways.add(tokens[0]);
            }
        }
        fu.close();
        return pathways;
    }
    
    private Set<GraphComponent> randomPickupComponents(String compFileName,
                                                       int idSize) throws IOException {
        List<GraphComponent> comps = loadComponentsInOrder(compFileName);
        Set<GraphComponent> picked = new HashSet<GraphComponent>();
        Set<String> ids = new HashSet<String>();
        while (ids.size() < idSize) {
            int index = (int) (Math.random() * comps.size());
            GraphComponent comp = comps.get(index);
            if (comp.getAllNodes().size() < CUT_OFF_OF_COMP_SIZE)
                continue;
            if (!picked.contains(comp)) {
                picked.add(comp);
                ids.addAll(comp.getAllNodes());
            }
        }
        return picked;
    }
    
    @Test
    public void randomPickupIdsViaFIsForCompTest() throws Exception {
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        if (totalGenes == null)
            loadTotalGenes();
        int permutation = 1000;
        List<String> geneList = new ArrayList<String>(totalGenes);
        int index = 0;
        for (int i = 0; i < permutation; i++) {
            index = (int) (Math.random() * geneList.size());
            String seed = geneList.get(index);
            // Random walk the seed to pick up components
            Set<String> sample = randomPickupComponent(seed, 
                                                      allFIs, 
                                                      20);
            GraphComponent comp = new GraphComponent();
            for (String gene : sample)
                comp.addNode(gene);
            double score = miCalculator.calculateScore(comp);
            System.out.println(score);
        }
    }
    
    private Set<String> randomPickupComponent(String seed,
                                              Set<String> fis,
                                              int size) {
        Set<String> comps = new HashSet<String>();
        Set<String> current = new HashSet<String>();
        current.add(seed);
        Set<String> next = new HashSet<String>();
        Set<String> checked = new HashSet<String>();
        while (comps.size() < size) {
            for (String tmp : current) {
                comps.add(tmp);
                if (comps.size() == size)
                    break;
                if (checked.contains(tmp))
                    continue;
                Set<String> interactors = getFIInteractors(tmp, fis);
                next.addAll(interactors);
                checked.add(tmp);
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
        return comps;
    }
    
    /**
     * Random pickup genes from the network and calculate their MI scores.
     * @throws Exception
     */
    @Test
    public void randomPickupIdsForCompTest() throws Exception {
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        if (totalGenes == null)
            loadTotalGenes();
        int permutation = 1000;
        List<String> geneList = new ArrayList<String>(totalGenes);
        for (int i = 0; i < permutation; i++) {
            Set<String> sample = sampleList(geneList,
                                            20);
            GraphComponent comp = new GraphComponent();
            for (String gene : sample)
                comp.addNode(gene);
            double score = miCalculator.calculateScore(comp);
            System.out.println(score);
        }
    }
    
    /**
     * This is a permutation test by shuffle ids to values mapping and do search component.
     * @throws Exception
     */
    @Test
    public void randomIdsSearchComponentsBasedOnMI() throws Exception {
        long time1 = System.currentTimeMillis();
        GraphComponentSearchEngine searcher = new GraphComponentSearchEngine();
        searcher.setMaxDepth(1);
        searcher.setSearchDepth(1);
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        String fileName = dirName + "Components1_1_GSE2023_ID_Permutated.txt";
        
//        // Search components based on van de Vijver data set
//        miCalculator.setSampleToPhenotypeInfo(dataAnalyzer.loadNejmPatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
//        String fileName = dirName + "Components2_1_Nejm_NO_P_I.txt";
        // Search based on GSE4922
//        miCalculator.setSampleToPhenotypeInfo(dataAnalyzer.loadGSE4922PatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesNormZScore.txt");
//        String fileName = dirName + "Components2_1_GSE4922_ID_Permutated.txt";  
        // Do randomization
        Map<String, List<Double>> idToValues = miCalculator.getIdToValues();
        Map<String, List<Double>> permutatedIdToValues = permutateIds(idToValues);
        miCalculator.setIdToValues(permutatedIdToValues);        
        searcher.setScoreCalculator(miCalculator);
        // Want to use the biggest component only
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + "FI73InGene_062008_NO_P_I_BigComp.txt";
        searcher.search(intFileName);
        List<GraphComponent> components = searcher.getFoundComponents();
        System.out.println("Total componenents: " + components.size());
        long time2 = System.currentTimeMillis();
        System.out.println("Time used: " + (time2 - time1));
        outputFoundComponents(fileName, components);
    }
    
    @Test
    public void randomSamplesSearchComponentsBasedOnMI() throws Exception {
        long time1 = System.currentTimeMillis();
        GraphComponentSearchEngine searcher = new GraphComponentSearchEngine();
        searcher.setMaxDepth(1);
        searcher.setSearchDepth(1);
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        String fileName = dirName + "Components1_1_GSE2023_Sample_Permutated.txt";
        
//        // Search components based on van de Vijver data set
//        miCalculator.setSampleToPhenotypeInfo(dataAnalyzer.loadNejmPatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
//        String fileName = dirName + "Components2_1_Nejm_NO_P_I.txt";
        // Search based on GSE4922
//        miCalculator.setSampleToPhenotypeInfo(dataAnalyzer.loadGSE4922PatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesNormZScore.txt");
//        String fileName = dirName + "Components2_1_GSE4922_ID_Permutated.txt";  
        // Do randomization
        List<String> samples = miCalculator.getSamples();
        List<String> permutatedSamples = permutateList(samples);
        miCalculator.setSamples(permutatedSamples);
        searcher.setScoreCalculator(miCalculator);
        // Want to use the biggest component only
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + "FI73InGene_062008_NO_P_I_BigComp.txt";
        searcher.search(intFileName);
        List<GraphComponent> components = searcher.getFoundComponents();
        System.out.println("Total componenents: " + components.size());
        long time2 = System.currentTimeMillis();
        System.out.println("Time used: " + (time2 - time1));
        outputFoundComponents(fileName, components);
    }
    
    private Map<String, List<Double>> permutateIds(Map<String, 
                                                   List<Double>> idToValues) {
        List<String> ids = new ArrayList<String>(idToValues.keySet());
        List<List<Double>> valuesList = new ArrayList<List<Double>>(idToValues.values());
        Map<String, List<Double>> permutated = new HashMap<String, List<Double>>();
        ObjectArrayList shuffleKeyList = new ObjectArrayList();
        shuffleKeyList.addAllOf(ids);
        shuffleKeyList.shuffle();
        ObjectArrayList shuffleValueList = new ObjectArrayList();
        shuffleValueList.addAllOf(valuesList);
        for (int i = 0; i < ids.size(); i++) {
            String id = (String) shuffleKeyList.get(i);
            List<Double> values = (List<Double>) shuffleValueList.get(i);
            permutated.put(id, values);
        }
        return permutated;
    }

    /**
     * This method is used to output a list of GraphComponent to an external file.
     * @param fileName
     * @param components
     * @throws IOException
     */
    public void outputFoundComponents(String fileName,
                                      List<GraphComponent> components) throws IOException {
        // Want to print out all components
        StringBuilder builder = new StringBuilder();
        fu.setOutput(fileName);
        int index = 1;
        for (GraphComponent comp : components) {
            builder.append(index).append("\t").append(comp.getScore());
            builder.append("\t").append(comp.getAllNodes().size());
            for (String node : comp.getAllNodes()) {
                builder.append("\t").append(node);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
            index ++;
        }
        fu.close();
    }
    
    @Test
    public void networkClusterTopMIGenes() throws Exception {
        String fileName = dirName + "GSE2034GeneToMI.txt";
        // Load all genes with mi >= 0.03
        List<String> genes = new ArrayList<String>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double mi = Double.parseDouble(tokens[1]);
            if (mi >= 0.03)
                genes.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total top genes: " + genes.size());
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterAnalyzer.cluster(genes, 0.1);
        int c = 0;
        // Calculate mutation information
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        for (Set<String> cluster : clusters) {
            GraphComponent graphComponent = new GraphComponent();
            for (String gene : cluster)
                graphComponent.addNode(gene);
            double mi = miCalculator.calculateScore(graphComponent);
            System.out.println(c + "\t" + cluster.size() + "\t" + mi);
            c ++;
        }
    }
    
    @Test
    public void sortGenesBasedOnMI() throws Exception {
        long time1 = System.currentTimeMillis();
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        Map<String, List<Double>> gene2Values = miCalculator.getIdToValues();
        final Map<String, Double> geneToMi = new HashMap<String, Double>();
        for (String gene : gene2Values.keySet()) {
            double mi = miCalculator.calculateScore(gene);
            geneToMi.put(gene, mi);
        }
        // Sorting based on mi
        List<String> geneList = new ArrayList<String>(geneToMi.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double mi1 = geneToMi.get(gene1);
                Double mi2 = geneToMi.get(gene2);
                return mi2.compareTo(mi1);
            }
        });
        // Print out the first 2000
        int c = 0;
        for (String gene : geneList) {
            System.out.println(gene + "\t" + geneToMi.get(gene));
            c ++;
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time used: " + (time2 - time1));
    }
    
    @Test
    public void generateEdgeWeightForFIsInFile() throws Exception {
        //Map<String, Double> fiToWeight = generateEdgeWeightForFIs();
        Map<String, Double> fiToWeight = calculateEdgeWeightFromTScore();
        // Need to output as a format that can be imported by MCL clustering
        //String outFileName = dirName + "FIsWithWeightGSE4922_091409.txt";
        String outFileName = dirName + "FIsWithWeightFromTTestGSE4922_091509.txt";
        fu.setOutput(outFileName);
        int unknown = 0;
        for (String fi : fiToWeight.keySet()) {
            Double weight = fiToWeight.get(fi);
            fu.printLine(fi + "\t" + weight);
            if (Math.abs(weight - 0.0) < 1.0E-6)
                unknown ++;
        }
        fu.close();
        System.out.println("Unknown: " + unknown + " out of " + fiToWeight.size());
    }
    
    private Map<String, Double> calculateEdgeWeightFromTScore() throws Exception {
        // Initialize an object for MI calculation
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
//        // Based on Wang et al data set
//        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
//        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
        // GSE4922
        miCalculator.setSampleToPhenotypeInfo(new BreastCancerArrayDataAnalyzer().loadGSE4922PatientInfo());
        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore_091409.txt");
        // Load FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        // Calculate MIs for each genes
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        for (String gene : genesInFIs) {
            Double mi = miCalculator.calculateTScore(gene);
            geneToScore.put(gene, mi);
        }
        Map<String, Double> fiToWeight = new HashMap<String, Double>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            Double score1 = geneToScore.get(gene1);
            String gene2 = fi.substring(index + 1);
            Double score2 = geneToScore.get(gene2);
            Double score = 0.0;
            if (score1.equals(Double.NEGATIVE_INFINITY) || score2.equals(Double.NEGATIVE_INFINITY)) {
                fiToWeight.put(fi, score);
                continue;
            }
            // Edge score from t-test is defined as average - sd.
            // If score less than 0, it will be round to 0.
            double avg = (score1 + score2) / 2.0;
            double sd = Math.sqrt(((score1 - avg) * (score1 - avg) + (score2 - avg) * (score2 - avg)) / 2.0);
            score = avg - sd;
            if (score < 0.0)
                score = 0.0;
            fiToWeight.put(fi, score);
        }
        return fiToWeight;
    }
    
    /**
     * This method is used to generate a random edge weight file for testing.
     * @throws Exception
     */
    @Test
    public void generateRandomEdgeWeightForFIsInFile() throws Exception {
        String fileName = dirName + "FIsWithWeightGSE4922_091409.txt";
        // Load all values 
        List<String> weights = new ArrayList<String>();
        fu.setInput(fileName);
        String line;
        List<String> fis = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            fis.add(tokens[0] + "\t" + tokens[1]);
            weights.add(tokens[2]);
        }
        fu.close();
        // Do a double permutations
        RandomDataImpl randomizer = new RandomDataImpl();
        Object[] randomWeights = randomizer.nextSample(weights, weights.size());
        Object[] randomFIs = randomizer.nextSample(fis, fis.size());
        String outFileName = dirName + "FIsWithRandomWeightsGSE4922_091409.txt";
        fu.setOutput(outFileName);
        for (int i = 0; i < randomWeights.length; i++) {
            fu.printLine(randomFIs[i] + "\t" + 
                         randomWeights[i]);
        }
        fu.close();
    }
    
    /**
     * This method is used to test some features in mutual information.
     * @throws Exception
     */
    @Test
    public void checkMIsForGenes() throws Exception {
        // Initialize an object for MI calculation
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
//        // Based on Wang et al data set
//        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
//        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
        // GSE4922
        miCalculator.setSampleToPhenotypeInfo(new BreastCancerArrayDataAnalyzer().loadGSE4922PatientInfo());
        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore_091409.txt");

        // Load FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        // Calculate MIs for each genes
        double max = 0.0d;
        for (String gene : genesInFIs) {
            Double mi = miCalculator.calculateScore(gene);
            if (mi > max)
                max = mi;
        }
        System.out.println("Max MI: " + max);
        // Calculate correlation for each fi
        Map<String, List<Double>> geneToValues = miCalculator.getIdToValues();
        Map<String, Double> fiToWeight = new HashMap<String, Double>();
        // Calculate Pearson correlation coefficient
        int index = 0;
        double maxCC = 0.0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            List<Double> values1 = geneToValues.get(gene1);
            List<Double> values2 = geneToValues.get(gene2);
            if (values1 == null || values2 == null) {
                continue;
            }
            Double cc = MathUtilities.calculatePearsonCorrelation(values1, values2);
            // Want to use the positive value only 
            if (cc < 0)
                cc = -cc;
            if (maxCC < cc)
                maxCC = cc;
        }
        System.out.println("Max CC: " + maxCC);
    }
    
    /**
     * This method is used to assign edge weights to FIs based on mutual information. An edge weight
     * is calculate as following:
     * 1). A FI has two genes. Calculate mutual information (MI) for each gene (mi1, mi2) based on
     * expression values and sample types
     * 2). Calculate Pearson correlation for these two genes (cc)
     * 3). The weight is defined as w = (mi1 + mi2 + cc) / 3.
     * 4). Note: Pearson in 2 probably should be changed to mi too as other two metrics. 
     * @throws Exception
     */
    private Map<String, Double> generateEdgeWeightForFIs() throws Exception {
        // Initialize an object for MI calculation
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
//        // Based on Wang et al data set
//        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
//        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
        // GSE4922
        miCalculator.setSampleToPhenotypeInfo(new BreastCancerArrayDataAnalyzer().loadGSE4922PatientInfo());
        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore_091409.txt");
        // Load FIs
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> genesInFIs = InteractionUtilities.grepIDsFromInteractions(fis);
        // Calculate MIs for each genes
        Map<String, Double> geneToMI = new HashMap<String, Double>();
        for (String gene : genesInFIs) {
            Double mi = miCalculator.calculateScore(gene);
            geneToMI.put(gene, mi);
        }
        Map<String, List<Double>> geneToValues = miCalculator.getIdToValues();
        Map<String, Double> fiToWeight = new HashMap<String, Double>();
        // This value is calculated from checkMIsForGenes()
        double maxMI = 0.0587199032343344;
        // Calculate Pearson correlation coefficient
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            List<Double> values1 = geneToValues.get(gene1);
            List<Double> values2 = geneToValues.get(gene2);
            Double cc = null;
            if (values1 == null || values2 == null) {
                cc = Double.NEGATIVE_INFINITY;
            }
            else {
                cc = MathUtilities.calculatePearsonCorrelation(values1, values2);
                // Want to use the positive value only 
                if (cc < 0)
                    cc = -cc;
            }
            Double score = ((geneToMI.get(gene1) + geneToMI.get(gene2)) / maxMI + cc) / 3.0;
            // Double score = (geneToMI.get(gene1) + geneToMI.get(gene2)) / 2.0d;
            if (score.equals(Double.NEGATIVE_INFINITY))
                score = 0.0d;
            fiToWeight.put(fi, score);
        }
        return fiToWeight;
    }
    
    /**
     * This method is used to search components in the FI network.
     * @throws Exception
     */
    @Test
    public void searchComponentsBasedOnMI() throws Exception {
        long time1 = System.currentTimeMillis();
        GraphComponentSearchEngine searcher = new GraphComponentSearchEngine();
        searcher.setMaxDepth(2);
        searcher.setSearchDepth(1);
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        //String fileName = dirName + "Components1_1_GSE2023.txt";
        String fileName = dirName + "Components2_1_GSE2034_062609.txt";
        
//        // Search components based on van de Vijver data set
//        miCalculator.setSampleToPhenotypeInfo(dataAnalyzer.loadNejmPatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
//        String fileName = dirName + "Components2_1_Nejm_NO_P_I.txt";
        // Search based on GSE4922
        // GSE4922
//        miCalculator.setSampleToPhenotypeInfo(loadGSE4922PatientInfo());
//        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesNormZScore.txt");
//        String fileName = dirName + "Components2_1_GSE4922_NO_P_I.txt";  
        
        searcher.setScoreCalculator(miCalculator);
        // Want to use the biggest component only
        //String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008_BigComp.txt";
        //String intFileName = R3Constants.RESULT_DIR + "FI73InGene_062008_NO_P_I_BigComp.txt";
        String intFileName = R3Constants.GENE_FI_BIG_COMP_FILE_NAME;
        searcher.search(intFileName);
        List<GraphComponent> components = searcher.getFoundComponents();
        System.out.println("Total componenents: " + components.size());
        long time2 = System.currentTimeMillis();
        System.out.println("Time used: " + (time2 - time1));
        outputFoundComponents(fileName, components);
    }
    
    /**
     * This method is used to calculate mutual information for clusters.
     * @throws IOException
     */
    @Test
    public void calculateClusterMutualInformation() throws Exception {
        long time1 = System.currentTimeMillis();
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredZScore.txt");
        // Try to use a randomized sample list
        List<String> originalSamples = miCalculator.getSamples();
        List<String> samples = permutateList(originalSamples);
        miCalculator.setSamples(samples);
        //Map<Integer, Set<String>> clusterToGenes = getClusterToGenes();
        //String fileName = dirName + "GSE2034ClusterMutInfo_New.txt";
        String fileName = dirName + "Components2_1_MI_random.txt";
        fu.setOutput(fileName);
        String compFileName = dirName + "Components2_1.txt";
        Map<Integer, Set<String>> clusterToGenes = getComponentToGenes(compFileName);
        for (Iterator<Integer> it = clusterToGenes.keySet().iterator(); it.hasNext();) {
            Integer cluster = it.next();
            Set<String> genes = clusterToGenes.get(cluster);
            if (genes.size() < 4)
                continue;
            GraphComponent component = new GraphComponent();
            for (String gene : genes)
                component.addNode(gene);
            double mi = miCalculator.calculateScore(component);
            fu.printLine(cluster + "\t" + mi);
        }
        fu.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
    }
    
    /**
     * This method is used to check scores in the clusterd generated from MCL.
     * @throws Exception
     */
    @Test
    public void calculateMCLClusterScore() throws Exception {
        long time1 = System.currentTimeMillis();
        Map<String, Double> fiToWeight = generateEdgeWeightForFIs();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        // To calcualte MIs
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
//        // Wang et al
//        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
//        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
        // GSE4922
        miCalculator.setSampleToPhenotypeInfo(new BreastCancerArrayDataAnalyzer().loadGSE4922PatientInfo());
        miCalculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore_091409.txt");

        // Try to use a randomized sample list
        String fileName = dirName + "FIsWithRandomWeightsMCLClustersGSE4922_I80_091409.txt";
        List<Set<String>> clusters = new MCLResultsAnalyzer().loadMCLClusters(fileName);
        int index = -1;
        System.out.println("Cluster\tScore\tMutual_Information\tSize\tFI_Size\tGenes");
        for (Set<String> cluster : clusters) {
            index ++;
            if (cluster.size() < 5)
                continue;
            Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
            if (fisInCluster.size() == 0) {
                System.out.println(index + ": " + cluster);
                continue;
            }
            double total = 0.0d;
            for (String fi : fisInCluster)
                total += fiToWeight.get(fi);
            double score = total / fisInCluster.size();
            List<String> geneList = new ArrayList<String>(cluster);
            Collections.sort(geneList);
            String genes = StringUtils.join(",", geneList);
            GraphComponent component = new GraphComponent();
            for (String gene : cluster)
                component.addNode(gene);
            double mi = miCalculator.calculateScore(component);
            System.out.println(index + "\t" + 
                               score + "\t" + 
                               mi + "\t" +
                               cluster.size() + "\t" + 
                               fisInCluster.size() + "\t" + 
                               genes);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
    }
    
    private List<GraphComponent> loadMCLClusterFromTScore(String fileName, 
                                                          MutualInformationScoreCalculator miCalculator,
                                                          double cutoff) throws Exception {
        List<GraphComponent> rtn = new ArrayList<GraphComponent>();
        List<Set<String>> clusters = new MCLResultsAnalyzer().loadMCLClusters(fileName);
        int index = -1;
        for (Set<String> cluster : clusters) {
            index ++;
            if (cluster.size() < 5)
                continue;
            double clusterScore = miCalculator.calculateTScore(cluster);
            if (clusterScore >= cutoff) {
                GraphComponent comp = new GraphComponent();
                comp.setId(index);
                for (String gene : cluster)
                    comp.addNode(gene);
                rtn.add(comp);
            }
        }
        return rtn;
    }
    
    @Test
    public void calculateMCLClusterTScore() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        String fileName = dirName + "FIsWithWeightFromTTestGSE2034_091509.txt";
        Map<String, Double> fiToWeight = loadFIsWithWeights(fileName);
        MCLResultsAnalyzer mclAnalyzer = new MCLResultsAnalyzer();
        List<Set<String>> clusters = mclAnalyzer.loadMCLClusters(dirName + "FIsWithWeightFromTTestMCLClustersGSE2034_I80_091509.txt");
        int index = 0;
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Wang et al
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore_091409.txt");
        System.out.println("Cluster\tCluster_Size\tCluster_FIs\tAvg_FI_Weight\tCluster_TScore\tHighest_TScore");
        for (Set<String> cluster : clusters) {
            if (cluster.size() < 5)
                continue;
            Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
            double total = 0.0;
            for (String fi : fisInCluster) {
                total += fiToWeight.get(fi);
            }
            double avg = total / fisInCluster.size();
            double clusterScore = miCalculator.calculateTScore(cluster);
            // Want to get the highest T_Score
            double highest = 0.0;
            for (String gene : cluster) {
                double tmp = miCalculator.calculateTScore(gene);
                if (tmp > highest)
                    highest = tmp;
            }
            System.out.println(index + "\t" + cluster.size() + "\t" +
                               fisInCluster.size() + "\t" + avg + "\t" + clusterScore + 
                               "\t" + highest);
            index ++;
        }
    }
    
}
