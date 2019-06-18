/*
 * Created on Apr 18, 2012
 *
 */
package org.reactome.cancer;

import static org.reactome.r3.util.R3Constants.BREAST_DIR;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;
/**
 * Some analysis methods related to paper published in Venet D, Dumont JE, Detours V (2011) 
 * Most Random Gene Expression Signatures Are Significantly Associated with Breast 
 * Cancer Outcome. PLoS Comput Biol 7(10): e1002240. doi:10.1371/journal.pcbi.1002240.
 * @author gwu
 *
 */
public class BreastCancerPCNAAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public BreastCancerPCNAAnalyzer() {
    }
    
    /**
     * Calculate a simple p-value score across five data sets. The p-value
     * score is the average of -logarithm.
     * @throws Exception
     */
    @Test
    public void calculatePValueScores() throws Exception {
        String fileName = BREAST_DIR + "PloSCompBioRandomSet/CoxPHForSignatures.txt";
        String outFileName = BREAST_DIR + "PloSCompBioRandomSet/CoxPHForSignaturesWithScores.txt";
        fu.setInput(fileName);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        fu.printLine(line);
        line = fu.readLine();
        fu.printLine(line + "\tP_Value_Score\tIs_Significant");
        int[] indices = new int[]{2, 4, 6, 8, 10};
        // Flag to check if there is any p-value > 0.05. If p-value > 0.05,
        // the total score will be 0.
        boolean isSignficant = false;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            isSignficant = true;
            double total = 0.0d;
            for (int index : indices) {
                Double pvalue = new Double(tokens[index]);
                if (pvalue > 0.05) {
                    isSignficant = false;
                }
                total += -Math.log10(pvalue);
            }
            fu.printLine(line + "\t" + total / indices.length + "\t" + isSignficant);
        }
        fu.close();
    }
    
    /**
     * This method is used to check relationships between module2 and other published gene
     * signatures, which have been collected by the paper.
     * @throws Exception
     */
    @Test
    public void checkModule2AndPublishedSignatures() throws Exception {
        Map<String, Set<String>> nameToSignatures = loadCancerSignatures();
//        for (String name : nameToSignatures.keySet()) {
//            Set<String> signature = nameToSignatures.get(name);
////            System.out.println(name + "\t" + signature.size()  + "\t" + signature);
//        }
        
        BreastMCLGeneExpModulePermutationTester breastAnalyzer = new BreastMCLGeneExpModulePermutationTester();
        String[] geneExpFileNames = breastAnalyzer.getSourceFileNames();
        
        Set<String> module2 = breastAnalyzer.getModule2();
        System.out.println("Gene overlapping:");
        System.out.println("Name\tSize\tShared\tP_Value");
        for (String name : nameToSignatures.keySet()) {
            Set<String> set = nameToSignatures.get(name);
            Set<String> shared = InteractionUtilities.getShared(module2, set);
            double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                                        set.size(), 
                                                                        module2.size(), 
                                                                        shared.size());
            System.out.println(name + "\t" +
                               set.size() + "\t" +
                               shared.size() + "\t" +
                               pvalue);
        }
        
        System.out.println("\nCheck gene expressionc correlation...");
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        List<Set<String>> list = new ArrayList<Set<String>>();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            System.out.println(geneExpFileNames[i]);
            System.out.println("Name\tSize\tCorrelation\tP_Value");
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(BREAST_DIR + geneExpFileNames[i]);
            for (String name : nameToSignatures.keySet()) {
                Set<String> set = nameToSignatures.get(name);
                list.clear();
                list.add(set);
                list.add(module2);
                PearsonsCorrelation correlation = calculateCorrelation(list, arrayHelper, geneToSampleToValue);
                System.out.println(name + "\t" + 
                                   set.size() + "\t" + 
                                   correlation.getCorrelationMatrix().getEntry(0, 1) + "\t" + 
                                   correlation.getCorrelationPValues().getEntry(0, 1));
            }
            System.out.println();
        }
    }
    
    /**
     * This method is used to calculate expression correlation between module2 and
     * meta-PCNA.
     * @throws Exception
     */
    @Test
    public void checkExpressionCorrelation() throws Exception {
        BreastMCLGeneExpModulePermutationTester breastAnalyzer = new BreastMCLGeneExpModulePermutationTester();
        String[] geneExpFileNames = breastAnalyzer.getSourceFileNames();
        
        List<Set<String>> list = new ArrayList<Set<String>>();
        list.add(loadPCNAGenes());
        list.add(breastAnalyzer.getModule2());
        
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            System.out.println(geneExpFileNames[i]);
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(BREAST_DIR + geneExpFileNames[i]);
            PearsonsCorrelation correlation = calculateCorrelation(list,
                                                                   arrayHelper,
                                                                   geneToSampleToValue);
            System.out.println("Correlation: " + correlation.getCorrelationMatrix().getEntry(0, 1));
            System.out.println("P-value: " + correlation.getCorrelationPValues().getEntry(0, 1));
            System.out.println();
        }
    }

    private PearsonsCorrelation calculateCorrelation(List<Set<String>> list,
                                                     CancerGeneExpressionCommon arrayHelper,
                                                     Map<String, Map<String, Double>> geneToSampleToValue) {
        List<Map<String, Double>> sampleToModuleScores = arrayHelper.generateSampleToGeneExpForClusters(geneToSampleToValue,
                                                                                                        list);
        List<Double> values1 = new ArrayList<Double>();
        List<Double> values2 = new ArrayList<Double>();
        Map<String, Double> map1 = sampleToModuleScores.get(0);
        Map<String, Double> map2 = sampleToModuleScores.get(1);
        for (String sample : map1.keySet()) {
            Double value1 = map1.get(sample);
            Double value2 = map2.get(sample);
            values1.add(value1);
            values2.add(value2);
        }
        PearsonsCorrelation correlation = MathUtilities.constructPearsonCorrelation(values1, values2);
        return correlation;
    }
    
    /**
     * This method is used to do survival analysis using the meta-PCNA genes.
     * @throws Exception
     */
    @Test
    public void doSurvivalAnalysisForMetaCNA() throws Exception {
        boolean includeModule2 = true;
        BreastMCLGeneExpModulePermutationTester breastAnalyzer = new BreastMCLGeneExpModulePermutationTester();
        String[] geneExpFileNames = breastAnalyzer.getSourceFileNames();
        String[] clinFileNames = breastAnalyzer.getClinFileNames();
        List<Map<String, Map<String, Double>>> geneToSampleToValues = new ArrayList<Map<String,Map<String,Double>>>();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(BREAST_DIR + geneExpFileNames[i]);
            geneToSampleToValues.add(geneToSampleToValue);
//            break;
        }
        List<Set<String>> list = new ArrayList<Set<String>>();
        Set<String> pcnaGenes = loadPCNAGenes();
        list.add(pcnaGenes);
        
        // For all signatures
        Map<String, Set<String>> nameToSignature = loadCancerSignatures();
        list.clear();
//        list.add(nameToSignature.get("YU")); // It has the highest p-value for expression correlation.
        list.add(nameToSignature.get("HALLSTROM")); // It has the medium p-value for expression correlation.
        
        if (includeModule2) {
            Set<String> module2 = breastAnalyzer.getModule2();
            list.add(module2);
        }
        for (int i = 0; i < geneToSampleToValues.size(); i++) {
            Map<String, Map<String, Double>> geneToSampleToValue = geneToSampleToValues.get(i);
            String tmpFileName = R3Constants.TEMP_DIR + "/TmpModuleToGeneExp.txt";
            arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                        tmpFileName,
                                                        list);
            File file = new File(tmpFileName);
            String clinFileName = BREAST_DIR  + clinFileNames[i];
            String[] output = breastAnalyzer.doSurvivalAnalysis(file,
                                                                new File(clinFileName), 
                                                                "coxph",
                                                                includeModule2 ? null : 0);
            System.out.println(geneExpFileNames[i]);
            System.out.println(output[0]);
            System.out.println(output[1]);
        }
        
//        if (true)
//            return;
        list.clear();
        int index = 1;
        List<String> sigNames = new ArrayList<String>();
        for (String name : nameToSignature.keySet()) {
            System.out.println(index + "\t" + name);
            list.add(nameToSignature.get(name));
            index ++;
            sigNames.add(name);
        }
        if (includeModule2) {
            list.add(breastAnalyzer.getModule2());
            sigNames.add("Module2");
        }
        String[] datasets = breastAnalyzer.getDataSetNames();
        for (int i = 0; i < geneToSampleToValues.size(); i++) {
            Map<String, Map<String, Double>> geneToSampleToValue = geneToSampleToValues.get(i);
            String tmpFileName = BREAST_DIR + "/PloSCompBioRandomSet/" + datasets[i] + "SampleToSignatures.txt";
//            String tmpFileName = R3Constants.TEMP_DIR + "/TmpModuleToGeneExp.txt";
            arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                        tmpFileName,
                                                        list,
                                                        sigNames);
            File file = new File(tmpFileName);
            String clinFileName = BREAST_DIR  + clinFileNames[i];
            String[] output = breastAnalyzer.doSurvivalAnalysis(file,
                                                                new File(clinFileName), 
                                                                "coxph",
                                                                null);
            System.out.println(geneExpFileNames[i]);
            System.out.println(output[0]);
            System.out.println(output[1]);
        }
    }
    
    /**
     * Use this method to load meta-PCNA genes.
     * @return
     * @throws IOException
     */
    private Set<String> loadPCNAGenes() throws IOException {
        String fileName = BREAST_DIR + "PloSCompBioRandomSet/super_PCNA_genes.txt";
        fu.setInput(fileName);
        String line = null;
        Set<String> genes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            String[] tokens = line.split(" ");
            for (String token : tokens)
                genes.add(token);
        }
        fu.close();
        return genes;
    }
    
    /**
     * Load cancer signatures tested by the paper. The original list was extracted from the R code.
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadCancerSignatures() throws IOException {
        Map<String, Set<String>> nameToSignatures = new HashMap<String, Set<String>>();
        String fileName = BREAST_DIR + "PloSCompBioRandomSet/cancer_signatures.txt";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null){
            if (line.startsWith("#"))
                continue;
            int index = line.indexOf(":");
            String name = line.substring(0, index);
            line = line.substring(index + 1).trim();
            String[] tokens = line.split(" ");
            Set<String> set = new HashSet<String>();
            for (String token : tokens)
                set.add(token.toUpperCase());
            nameToSignatures.put(name, set);
        }
        fu.close();
        return nameToSignatures;
    }
    
    /**
     * Check overlapping between meta-PCNA genes and our module 2.
     * @throws IOException
     */
    @Test
    public void checkOverlapping() throws Exception {
        BreastMCLGeneExpModulePermutationTester breastAnalyzer = new BreastMCLGeneExpModulePermutationTester();
        Set<String> module2 = breastAnalyzer.getModule2();
        System.out.println("Module 2 genes: " + module2.size());
        Set<String> pcnaGenes = loadPCNAGenes();
        System.out.println("PCNA genes: " + pcnaGenes.size());
        Set<String> shared = InteractionUtilities.getShared(module2, pcnaGenes);
        System.out.println("Shared genes: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES, 
                                                                    module2.size(),
                                                                    pcnaGenes.size(),
                                                                    shared.size());
        System.out.println("P-value: " + pvalue);
        System.out.println("Shared genes: " + shared);
    }
    
}
