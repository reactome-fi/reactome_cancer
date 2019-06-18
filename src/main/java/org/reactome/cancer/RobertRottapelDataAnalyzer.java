/*
 * Created on Jun 7, 2012
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to analyze data sets from Robert Rottapel's group.
 * @author gwu
 *
 */
public class RobertRottapelDataAnalyzer {
    private FileUtility fu = new FileUtility();
    private final String OVARIAN_DIR_NAME = CancerResequenceDataSetAnalyzer.OVARIAN_DIR_NAME;
    private final String DIR_NAME = OVARIAN_DIR_NAME + "Rottapel/";
    
    public RobertRottapelDataAnalyzer() {
    }
    
    /**
     * Generate a list of genes from MCL network clustering.
     * @throws Exception
     */
    @Test
    public void generateClusterGeneList() throws Exception {
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        String clusterFileName = DIR_NAME + "MCLClusterResults_061512.txt";
        List<Set<String>> clusters = clusterAnalyzer.loadNetworkClusters(clusterFileName);
        
        List<String> genes = new ArrayList<String>();
        for (Set<String> cluster : clusters) {
            genes.addAll(cluster);
        }
        Collections.sort(genes);
        
        for (String gene : genes) {
            System.out.println(gene);
        }
    }
    
    /**
     * MCL clustering analysis using p-values from the original file.
     * @throws Exception
     */
    @Test
    public void mclClusterAnalysis() throws Exception {
        String srcFileName = DIR_NAME + "ovarian_screens_logIC50_gene_panel_complete.txt";
        int pvalueIndex = 4;
//        String srcFileName = DIR_NAME + "ovarian_screens_SR_gene_panel_complete.txt";
//        int pvalueIndex = 6;
        Map<String, Double> geneToPValue = new HashMap<String, Double>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[1];
            if (gene.length() == 0)
                continue;
            Double value = new Double(tokens[pvalueIndex]); // p-values
            // Duplicated values for genes exist in the file. Pickup the best p-value
            Double oldValue = geneToPValue.get(gene);
            if (oldValue == null || oldValue > value)
                geneToPValue.put(gene, value);
        }
        fu.close();
        
        double valueCutoff = 1.3d; // two genes should have p-value <= 0.05
        //        double valueCutoff = 1.0d;
        int sizeCutoff = 3;
        
        MCLClusterWrapper clusterWrapper = new MCLClusterWrapper();
        clusterWrapper.setKeepTempFile(false);
        //        clusterWrapper.setInflation(4.5d);
        List<Set<String>> clusters = clusterWrapper.mclClusterForGeneScores(geneToPValue,
                                                                            true,
                                                                            sizeCutoff,
                                                                            valueCutoff);
        
        System.out.println("Total clusters: " + clusters.size());
        String clusterFileName = DIR_NAME + "MCLClusterResults_061512.txt";
        new NetworkClusterAnalyzer().outputNetworkClusters(clusters, clusterFileName);
//        List<Double> clusterValues = clusterWrapper.calculateAverageCorrelations(clusters,
//                                                                                 fis, 
//                                                                                 fisWithValues);
//        for (int i = 0; i < clusters.size(); i++) {
//            System.out.println(i + "\t" + clusters.get(i).size() + "\t" + 
//                               clusterValues.get(i) + "\t" + 
//                               clusters.get(i));
//        }
        
        System.out.println("\nSurvival analysis using gene expression data sets:");
        checkWithOVByGeneExp(clusters);
        
        System.out.println("\nSurvival analysis using TCGA mutation data:");
        Set<String> allClusterGenes = new HashSet<String>();
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            System.out.println("Cluster " + i);
            checkWithTCGAOVMutation(cluster);
            allClusterGenes.addAll(cluster);
        }
        System.out.println("All clustered genes merged together: ");
        checkWithTCGAOVMutation(allClusterGenes);
    }
    
    private Set<String> getGeneSignature(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Set<String> genes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[1].trim();
            if (gene.length() == 0 || gene.equals("-"))
                continue;
            genes.add(gene);
        }
        fu.close();
        return genes;
    }
    
    private void filterRibosomeGenes(Set<String> genes) {
        for (Iterator<String> it = genes.iterator(); it.hasNext();) {
            String gene = it.next();
            if (gene.startsWith("RPL"))
                it.remove();
        }
        // The following genes are related to mRNA synthesis using polymerase
        String[] extraGenes = new String[] {
                "POLR2E",
                "SF3B5",
                "U2AF2",
                "NHP2L1",
                "POLR2K",
                "PABPN1"
        };
        for (String gene : extraGenes)
            genes.remove(gene);
    }
    
    /**
     * Check the provided gene signautre against TCGA or other OV gene expression data sets
     * for survival analysis.
     * @throws Exception
     */
    @Test
    public void checkWithTCGADatasetsByGeneExpression() throws Exception {
        String fileName = DIR_NAME + "essential gene signature.txt";
        Set<String> rottapelGenes = getGeneSignature(fileName);
        System.out.println("Total genes in " + fileName + ": " + rottapelGenes.size());
        // Filter out ribosome genes
//        filterRibosomeGenes(rottapelGenes);
        System.out.println("\tAfter filtering out ribosome genes: " + rottapelGenes.size());
        // In order to use the following statement
        List<Set<String>> geneList = new ArrayList<Set<String>>();
        geneList.add(rottapelGenes);
        
        checkWithOVByGeneExp(geneList);
    }

    private void checkWithOVByGeneExp(List<Set<String>> geneList) throws IOException {
        OvarianMCLGeneExpModuleAnalyzer moduleAnalyzer = new OvarianMCLGeneExpModuleAnalyzer();
        String[] srcFileNames = moduleAnalyzer.getSourceFileNames();
        String[] clinFileNames = moduleAnalyzer.getClinFileNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        
        BreastMCLGeneExpModulePermutationTester survivalHelper = new BreastMCLGeneExpModulePermutationTester();
        
        //String tmpFileName = "tmp/TmpModuleToGeneExp.txt";
        for (int i = 0; i < srcFileNames.length; i++) {
            String expFileName = srcFileNames[i];
            File expFile = new File(expFileName);
            String clinFileName = clinFileNames[i];
            System.out.println("Source file: " + expFileName);
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
            String tmpFileName = DIR_NAME + "shRNA " + expFile.getName();
            arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                        tmpFileName,
                                                        geneList);
            File file = new File(tmpFileName);
            String[] output = survivalHelper.doSurvivalAnalysis(file,
                                                                new File(clinFileName), 
                                                                "coxph",
                                                                null);
            System.out.println(output[0]);
        }
    }
    
    /**
     * Check the provided gene signature against the TCGA OV data set for survival assuming these genes
     * are mutated genes.
     * @throws Exception
     */
    @Test
    public void checkWithTCGADatasetByMutation() throws Exception {
        String fileName = DIR_NAME + "essential gene signature.txt";
        Set<String> rottapelGenes = getGeneSignature(fileName);
        
        checkWithTCGAOVMutation(rottapelGenes);
    }

    private void checkWithTCGAOVMutation(Set<String> rottapelGenes) throws IOException {
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        
        TCGAOvarianCancerAnalyzer tcgaAnalyzer = new TCGAOvarianCancerAnalyzer();
        String clinFileName = OVARIAN_DIR_NAME + "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt";
        Map<String, Set<String>> sampleToGenes = tcgaAnalyzer.loadSampleToNonSynonymousMutatedGenes(false);

//        filterRibosomeGenes(rottapelGenes);
        
        tcgaAnalyzer.doSurvivalAnalysisForMutationModule(survivalHelper, 
                                                         sampleToGenes,
                                                         rottapelGenes, 
                                                         clinFileName,
                                                         "Rottapel");
    }
    
    /**
     * Compare to several TCGA gene signatures.
     * @throws Exception
     */
    @Test
    public void compareToTCGAGenes() throws Exception {
        String fileName = DIR_NAME + "essential gene signature.txt";
        Set<String> rottapelGenes = getGeneSignature(fileName);
        System.out.println("Total genes in " + fileName + ": " + rottapelGenes.size());
        fu.saveCollection(rottapelGenes, DIR_NAME + "essential signature genes.txt");
//        if (true)
//            return;
        // Get signature genes from TCGA
        Set<String> tcgaModuleGenes = fu.loadInteractions(OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt");
        System.out.println("TCGA Module 6 and 7 genes: " + tcgaModuleGenes.size());
        Set<String> shared = InteractionUtilities.getShared(rottapelGenes, tcgaModuleGenes);
        System.out.println("Shared genes: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES, 
                                                                    rottapelGenes.size(),
                                                                    tcgaModuleGenes.size(), 
                                                                    shared.size());
        System.out.println("p-value: " + pvalue);
        OvarianMCLGeneExpModuleAnalyzer mclAanalyzer = new OvarianMCLGeneExpModuleAnalyzer();
        Set<String> ovMCLGenes = mclAanalyzer.getSelectedGenesFromSuperpc();
        System.out.println("Total MCL genes: " + ovMCLGenes.size());
        shared = InteractionUtilities.getShared(rottapelGenes, ovMCLGenes);
        System.out.println("Shared genes: " + shared.size() + " " + shared);
        pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                             rottapelGenes.size(),
                                                             ovMCLGenes.size(),
                                                             shared.size());
        System.out.println("p-value: " + pvalue);
    }
    
    
}
