/*
 * Created on Jan 11, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.EntrezGeneAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class contains methods related to MCL gene expression module based analysis on the breast cancer
 * data set.
 * @author wgm
 *
 */
public class BreastMCLGeneExpModuleAnalyzer {
    protected final String DIR_NAME = R3Constants.BREAST_DIR;
    private final FileUtility fu = new FileUtility();
    
    public BreastMCLGeneExpModuleAnalyzer() {
    }
    
    /**
     * This method is used to load 76 gene signature from the GSE2034 data set
     * (Wang et al, 2005).
     * @return
     * @throws IOException
     */
    private Set<String> loadGSE2034GeneSignature() throws IOException {
        String fileName = DIR_NAME + "GSE2034/Signature_76_Genes.txt";
        fu.setInput(fileName);
        String line = null;
        Set<String> probesets = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            String[] tokens = line.split(" ");
            probesets.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total probesets: " + probesets.size());
        GSESoftTextDataHandler softHandler = new GSESoftTextDataHandler();
        Map<String, String> probesetToGene = softHandler.loadProbeIdToGene(DIR_NAME + "GSE2034/GSE2034_family.soft");
        Set<String> signature = new HashSet<String>();
        for (String probeset : probesets) {
            String gene = probesetToGene.get(probeset);
            if (gene != null) {
                String[] tokens = gene.split("(-|///)");
                for (String token : tokens)
                    signature.add(token.trim());
            }
        }
        return signature;
    }
    
    @Test
    public void testLoadGSE2034GeneSignature() throws Exception {
        Set<String> signature = loadGSE2034GeneSignature();
        System.out.println("Total genes in signature: " + signature.size());
        for (String gene : signature)
            System.out.println(gene);
    }
    
    /**
     * This method is used to check gene signature overlappings.
     * @throws IOException
     */
    @Test
    public void checkGeneOverlappingBetweenMCLModule2AndNejmSignatureGenes() throws Exception {
        String nejmSignatureFile = "datasets/BreastCancer/nejm_table2/Gene70_Partial.txt";
        Set<String> nejmSignature = fu.loadInteractions(nejmSignatureFile);
        // Signature from Wang et al (GSE2034)
        Set<String> signature = loadGSE2034GeneSignature();
        System.out.println("Total signature genes: " + signature.size());
        Set<String> sharedSignature = InteractionUtilities.getShared(signature, nejmSignature);
        System.out.println("Shared sigature: " + sharedSignature.size() + " " + sharedSignature);
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        signature = entrezAnalyzer.normalizeGeneNames(signature);
        System.out.println("After name normalization: " + signature.size());
//        for (String gene : signature)
//            System.out.println(gene);
        // Load MCL genes
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        List<Set<String>> clusters =  common.selectGeneExpClusters(clusterFileName,
                                                                   fiToCorFileName, 
                                                                   OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                                                   OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF);
//        Set<String> cluster2 = clusters.get(1);
//        System.out.println("Genes in cluster2: " + cluster2.size());
//        Set<String> shared = InteractionUtilities.getShared(signature, cluster2);
//        System.out.println("Shared: " + shared.size());
//        System.out.println("Shared genes: " + shared);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Module\tModule_Size\tShared\tP-value");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            Set<String> shared = InteractionUtilities.getShared(signature, cluster);
            double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(), 
                                                                        signature.size(),
                                                                        cluster.size(),
                                                                        shared.size());
            System.out.println((i + 1) + "\t" + cluster.size() + "\t" + 
                               shared.size() + "\t" + pvalue + "\t" + shared);
        }
    }

    
    @Test
    public void checkClinicalInformationForGSE3143() throws IOException {
        // Get map from sample names to GSM id
        Map<String, String> sampleNameToGSMId = new HashMap<String, String>();
        String srcFileName = DIR_NAME + "GSE3143/GSE3143_series_matrix.txt";
        fu.setInput(srcFileName);
        String line = null;
        List<String> sampleIds = new ArrayList<String>();
        List<String> gsmIds = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("!Sample_title")) {
                String[] tokens = line.split("\t");
                for (int i = 1; i < tokens.length; i++) {
                    String sampleId = tokens[i];
                    sampleId = sampleId.substring(1, sampleId.length() - 1);
                    sampleIds.add(sampleId);
                }
            }
            else if (line.startsWith("\"ID_REF\"")) {
                String[] tokens = line.split("\t");
                for (int i = 1; i < tokens.length; i++) {
                    String gsmId = tokens[i];
                    gsmId = InteractionUtilities.removeQuotationMarks(gsmId);
                    gsmIds.add(gsmId);
                }
            }
        }
        fu.close();
        for (int i = 0; i < sampleIds.size(); i++) {
            String sampleId = sampleIds.get(i);
            // Need some manipulation
            int index = sampleId.indexOf("-");
            String number = sampleId.substring(index + 1);
            if (sampleId.startsWith("KF")) {
                // Need three digits
                while (number.length() < 3) {
                    number = "0" + number;
                }
            }
            else if (sampleId.startsWith("T")) {
                // Need four digist
                while (number.length() < 4)
                    number = "0" + number;
            }
            String gsmId = gsmIds.get(i);
            sampleNameToGSMId.put(sampleId.substring(0, index) + "-" + number, 
                                  gsmId);
        }
        Collections.sort(sampleIds);
        for (String sampleId : sampleIds)
            System.out.println(sampleId);
        // Load clinical information data
        srcFileName = DIR_NAME + "GSE3143/Breast_Clinical.txt";
        fu.setInput(srcFileName);
        fu.setOutput(DIR_NAME + "GSE3143/GSE3143_Survival_Info.txt");
        fu.printLine("Sample\tERLev\tOSEVENT\tOSDURATION");
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gsmId = sampleNameToGSMId.get(tokens[0]);
            fu.printLine(gsmId + "\t" + 
                         tokens[1] + "\t" + 
                         tokens[2] + "\t" +
                         tokens[3]);
        }
        fu.close();
    }
    
    @Test
    public void generateFIsWithGeneExpCorr() throws IOException {
//        String srcFileName = DIR_NAME + "NejmLogRatioNormZScore.txt";
//        String outFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        String srcFileName = DIR_NAME + "NejmLogRatioNormGlobalZScore_070111.txt";
        String outFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_070111.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        helper.generateFIsWithGeneExpCorr(srcFileName,
                                          outFileName,
                                          null);
    }
    
    @Test
    public void generateNejmClinInfoFile() throws IOException {
        String srcFileName = DIR_NAME + "Nejm_ClinicalData_Table.txt";
        // Create a simplified version
        String targetFileName = DIR_NAME + "Nejm_Clin_Simple.txt";
        fu.setInput(srcFileName);
        fu.setOutput(targetFileName);
        // Note: the unit of survival times is year actually.
        fu.printLine("SampleId\tOSEVENT\tOSDURATION");
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sampleId = tokens[0];
            String osEvent = tokens[4];
            String osDurartion = tokens[5];
            fu.printLine("Sample " + sampleId + "\t" +
                         osEvent + "\t" + 
                         osDurartion);
        }
        fu.close();
    }
    
    @Test
    public void generateGSE4922ClinInfoFile() throws IOException {
        String srcFileName = DIR_NAME + "GSE4922/Patient_Info_Tab_Survival.txt";
        // Create a simplified version
//        String targetFileName = DIR_NAME + "GSE4922_Clin_Simple.txt";
        String targetFileName = DIR_NAME + "GSE4922_Clin_Full.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        fu.setOutput(targetFileName);
        // Note: the unit of survival times is year actually.
        // Note: OSEVENT Maybe re-recurrence
        //fu.printLine("SampleId\tOSEVENT\tOSDURATION");
        String line = fu.readLine();
        fu.printLine(line);
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sampleId = tokens[0];
            int index = sampleId.indexOf("/");
            sampleId = sampleId.substring(0, index);
//            String osEvent = tokens[8];
//            String osDurartion = tokens[7];
//            fu.printLine(sampleId + "\t" +
//                         osEvent + "\t" + 
//                         osDurartion);
            builder.append(sampleId);
            for (int i = 0; i < tokens.length; i++) {
                if (i == 0)
                    continue;
                String t = tokens[i];
                builder.append("\t");
                if (i == 12) {
                    if (t.equals("ER-"))
                        builder.append("0");
                    else if (t.equals("ER+"))
                        builder.append("1");
                    else
                        builder.append("");
                }
                else if (i == 13) {
                    if (t.equals("LN-"))
                        builder.append("0");
                    else if (t.equals("LN+"))
                        builder.append("1");
                    else
                        builder.append("");
                }
                else if (i == 14) {
                    if (t.equals("p53-"))
                        builder.append("0");
                    else if (t.equals("p53+"))
                        builder.append("1");
                    else
                        builder.append("");
                }
                else
                    builder.append(tokens[i]);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void generateGSE4922SampleToGeneExpClusters() throws IOException {
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        String output = DIR_NAME + "GSE4922SampleToGeneExpFIsWithNejmGeneExpAbsCorr_011111.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        String geneExpFileName = DIR_NAME + "GSE4922FilteredOnSamplesZScore_091409.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);
    }
    
    @Test
    public void generateGSE3143SampleToGeneExpClusters() throws IOException {
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
//        String output = DIR_NAME + "GSE3143/GSE3143SampleToGeneExpFIsWithNejmGeneExpAbsCorr_012111.txt";
//        String output = DIR_NAME + "GSE3143/GSE3143_log_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_012111.txt";
        String output = DIR_NAME + "GSE3143/GSE3143_z_on_samples_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_012111.txt";
//        String output = DIR_NAME + "GSE3143/GSE3143_log_z_on_samples_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_012111.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
//        String geneExpFileName = DIR_NAME + "GSE3143/GSE3143_MappedGenes_z.txt";
//        String geneExpFileName = DIR_NAME + "GSE3143/GSE3143_MappedGenes_log_z.txt";
        String geneExpFileName = DIR_NAME + "GSE3143/GSE3143_MappedGenes_z_on_samples_012111.txt";
//        String geneExpFileName = DIR_NAME + "GSE3143/GSE3143_MappedGenes_log_z_on_samples_012111.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);
    }
    
    @Test
    public void generateGSE1992SampleToGeneExpClusters() throws IOException {
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
//        String output = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_020711.txt";
//        String output = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_020711.txt";
//        String output = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_POE_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_042611.txt";
//        String output = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt";
//        String output = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt";
        String output = DIR_NAME + "GSE1992/GSE1992_GPL1390_Gene_Exp_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080411.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
//        String geneExpFileName = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_z_020711.txt";
//        String geneExpFileName = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_020711.txt";
//        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_020711_POE.txt";
//        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_z_080311.txt";
//        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_080311.txt";
        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_GPL1390_Gene_Exp_080411.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);

    }
    
    @Test
    public void generateGSE18229SampleToGeneExpClusters() throws IOException {
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        String outputs[] = new String[] {
//                "GSE18229_GPL1390_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt",
//                "GSE18229_GPL1390_Gene_Exp_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt",
//                "GSE18229_Gene_Exp_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt",
//                "GSE18229_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt"
                "GSE18229_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_111811.txt"
        };
        String geneExpFileNames[] = new String[] {
//                "GSE18229_GPL1390_Gene_Exp_z_080311.txt",
//                "GSE18229_GPL1390_Gene_Exp_080311.txt",
//                "GSE18229_Gene_Exp_080311.txt",
//                "GSE18229_Gene_Exp_z_080311.txt"
                "GSE18229_Gene_Exp_z_111811.txt"
        };
//        String output = DIR_NAME + "GSE18229/GSE18229_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt";
//        String output = DIR_NAME + "GSE18229/GSE18229_GPL1390_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_080311.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
//        String geneExpFileName = DIR_NAME + "GSE18229/GSE18229_Gene_Exp_z_080311.txt";
//        String geneExpFileName = DIR_NAME + "GSE18229/GSE18229_GPL1390_Gene_Exp_z_080311.txt";
        for (int i = 0; i < geneExpFileNames.length; i++) {
            String geneExpFileName = DIR_NAME + "GSE18229/" + geneExpFileNames[i];
            String output = DIR_NAME + "GSE18229/" + outputs[i];
            Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
            common.generateSampleToGeneExpClusters(clusterFileName,
                                                   fiToCorFileName, 
                                                   geneToSampleToValue, 
                                                   OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                                   OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                                   output);
        }
    }
    
    @Test
    public void generateGSE1992SampleToGeneExpClustersWithPlugInOutput() throws IOException {
        String clusterFileName = DIR_NAME + "MCLModulesOnNejmPOE_042611.txt";
        String output = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_POE_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_MCLModulesOnNejmPOE_042611.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
//        String geneExpFileName = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_z_020711.txt";
//        String geneExpFileName = DIR_NAME + "GSE1922/GSE1992_Gene_Exp_020711.txt";
        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_020711_POE.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        List<Set<String>> clusters = common.loadGeneExpClusters(clusterFileName);
        common.generateSampleToGeneExpClusters(geneToSampleToValue, 
                                               output,
                                               clusters);
    }
    
    @Test
    public void generateGSE1456SampleToGeneExpClusters() throws IOException {
        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        String output = DIR_NAME + "GSE1456/GSE1456_Gene_Exp_z_SampleToGeneExpFIsWithNejmGeneExpAbsCorr_020911.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        String geneExpFileName = DIR_NAME + "GSE1456/GSE1456_z_020911.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);

    }
    
    @Test
    public void generateNejmSampleToGeneExpClusters() throws IOException {
//        String clusterFileName = DIR_NAME + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
//        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
//        String output = DIR_NAME + "SampleToGeneExpFIsWithNejmGeneExpAbsCorr_011111.txt";
//        String geneExpFileName = DIR_NAME + "NejmLogRatioNormZScore.txt";
        
        
        String clusterFileName = DIR_NAME + "MCLModulesNejmLogRatioNormGlobalZScore_070111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_070111.txt";
        String output = DIR_NAME + "SampleToGeneExpFIsWithNejmGeneExpAbsCorr_070111.txt";
        String geneExpFileName = DIR_NAME + "NejmLogRatioNormGlobalZScore_070111.txt";

        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);
    }
    
    /**
     * This method is used to check overlapping between Module 2 genes (index 1) and 
     * 16 kinase from PMID 18245477.
     * @throws Exception
     */
    @Test
    public void checkOverlappingBetweenModule2AndPMID18245477Genes() throws Exception {
        Set<String> module = loadMCLModule2();
        System.out.println("Total genes in module: " + module.size());
//        // Check gene overlapping between 16 genes in PMID: 18245477
//        String geneNames = "AURKA AURKB BUB1 BUB1B CDC2 CDC7 CHEK1 MASTL " +
//                "MELK NEK2 PBK PLK1 PLK4 SRPK1 TTK VRK1";
//        String[] names = geneNames.split(" ");
//        Set<String> otherGenes = new HashSet<String>();
//        for (String name : names)
//            otherGenes.add(name);
        
        Set<String> otherGenes = loadOncotypeDXGenes(false);
        
//        String test = "datasets/BreastCancer/nejm_table2/Gene70.txt";
//        fu.setInput(test);
//        String line = fu.readLine();
//        fu.close();
//        String[] tokens = line.split("\t");
//        for (String token : tokens)
//            System.out.println(token);
//        
//        if (true)
//            return;
//        
//        String nejmSignatureFile = "datasets/BreastCancer/nejm_table2/Gene70_Partial.txt";
//        Set<String> otherGenes = fu.loadInteractions(nejmSignatureFile);
        System.out.println("Total genes: " + otherGenes.size());
        
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        otherGenes = entrezAnalyzer.normalizeGeneNames(otherGenes);
        System.out.println("After name normalization: " + otherGenes.size());
        
        Set<String> shared = InteractionUtilities.getShared(otherGenes, module);
        System.out.println("Shared: " + shared.size());
        // Try to get total genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        otherGenes.retainAll(totalGenes);
        System.out.println("in FIs: " + otherGenes.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(),
                                                                    otherGenes.size(),
                                                                    module.size(), 
                                                                    shared.size());
        System.out.println("pvalue from hyper-geometric: " + pvalue);
        System.out.println("Shared genes are: " + shared);
    }
    
    protected Set<String> loadOncotypeDXGenes(boolean needReferences) throws IOException {
        String fileName = DIR_NAME + "OncotypeDXGenes.txt";
        fu.setInput(fileName);
        String line = null;
        Set<String> genes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (!needReferences && tokens.length > 2 && tokens[2].equals("Reference"))
                continue;
            genes.add(tokens[1]);
        }
        fu.close();
        return genes;
    }

    /**
     * Load the second module which is most interested in the breast cancer analysis.
     * @return
     * @throws IOException
     */
    protected Set<String> loadMCLModule2() throws IOException {
        String fileName = DIR_NAME + "SelectedMCLModules_7_025.txt";
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        Set<String> module = common.loadMCLModule(fileName, 1);
        return module;
    }
    
    /**
     * This method is used to compare two sets of MCL modules: one without global z-score
     * transformation, another with z-score transformation.
     * @throws Exception
     */
    @Test
    public void compareTwoMCLModules() throws Exception {
        //      String clusterFileName = R3Constants.BREAST_DIR + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String clusterFileName = R3Constants.BREAST_DIR + "MCLModulesNejmLogRatioNormGlobalZScore_070111.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clustersNew = clusterHelper.loadNetworkClusters(clusterFileName);
        clusterFileName = R3Constants.BREAST_DIR + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String fiToCorFileName = DIR_NAME + "FIsWithNejmGeneExpAbsCorr_011111.txt";
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        List<Set<String>> clustersOld = common.selectGeneExpClusters(clusterFileName,
                                                                     fiToCorFileName, 
                                                                     OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF, 
                                                                     OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF);
        // Want to do a sorting
        Comparator<Set<String>> sorter = new Comparator<Set<String>>() {
            public int compare(Set<String> set1, Set<String> set2) {
                if (set1.size() > set2.size())
                    return -1;
                if (set1.size() < set2.size())
                    return 1;
                String text1 = set1.toString();
                String text2 = set2.toString();
                return text1.compareTo(text2);
            }
        };
        Collections.sort(clustersNew, sorter);
        Collections.sort(clustersOld, sorter);
        for (int i = 0; i < clustersOld.size(); i++) {
            Set<String> clusterNew = clustersNew.get(i);
            Set<String> clusterOld = clustersOld.get(i);
            boolean isSame = clusterNew.equals(clusterOld);
            System.out.println(i + ": " + isSame);
            if (!isSame) {
                System.out.println("Old: " + clusterOld);
                System.out.println("New: " + clusterNew);
            }
        }
    }
    
    /**
     * This method is used to generate samples to gene expression modules for genes
     * in module 2, which is significantly related to breast cancer patient overall
     * survival times.
     * @throws IOException
     */
    @Test
    public void generateGeneExpressionProfileForModule2() throws Exception {
        Set<String> module = loadMCLModule2();
        System.out.println("Total genes in module: " + module.size());
        // Nejm data set
//        String geneExpFileName = DIR_NAME + "NejmLogRatioNormZScore.txt";
//        String outFileName = DIR_NAME + "NejmLogRatioNormZScore_Module2_021511.txt";
        
        // GSE4922
//        String geneExpFileName = DIR_NAME + "GSE4922FilteredOnSamplesZScore_091409.txt";
//        String outFileName = DIR_NAME + "GSE4922FilteredOnSamplesZScore_091409_Module2_021511.txt";
        
        // GSE3143
//        String geneExpFileName = DIR_NAME + "GSE3143/GSE3143_MappedGenes_z_on_samples_012111.txt";
//        String outFileName = DIR_NAME + "GSE3143_MappedGenes_z_on_samples_012111_Module2_021511.txt";
        
        // GSE1992
//        String geneExpFileName = DIR_NAME + "GSE1992/GSE1992_Gene_Exp_020711.txt";
//        String outFileName = DIR_NAME + "GSE1992_Gene_Exp_020711_Module2_021511.txt";
        
        //GSE1456
//        String geneExpFileName = DIR_NAME + "GSE1456/GSE1456_z_020911.txt";
//        String outFileName = DIR_NAME + "GSE1456_z_020911_Module2_021511.txt";
//        
//        Map<String, Map<String, Double>> geneToSampleToValue = common.loadGeneExp(geneExpFileName);
//        Set<String> keySet = geneToSampleToValue.keySet();
//        for (Iterator<String> it = keySet.iterator(); it.hasNext();) {
//            String gene = it.next();
//            if (module.contains(gene))
//                continue;
//            it.remove();
//        }
//        common.outputGeneExp(geneToSampleToValue, outFileName);
    }
    
    /**
     * This method is used to generate a nice-formatted file for a list of luminal A cancers
     * based on a paper from pubmed 18245477. The source file was copy/paste from Supp table 3, 
     * which is a PDF file itself.
     * @throws IOException
     */
    @Test
    public void processLuminalAFromPubmed18245477() throws IOException {
        String srcFileName = DIR_NAME + "Nejm_LuminalA_on_pubmed_18245477.txt";
        fu.setInput(srcFileName);
        String line = null;
        Map<Integer, String> sampleToType = new HashMap<Integer, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("(\t)|( )");
            for (int i = 0; i < tokens.length; i++) {
                String token = tokens[i];
                if (token.matches("(\\d)+")) {
                    sampleToType.put(new Integer(token), tokens[++i]);
                }
            }
        }
        fu.close();
        System.out.println("Total samples: " + sampleToType.size());
        String targetFileName = DIR_NAME + "Nejm_LuminalA_on_pubmed_18245477_formatted.txt";
        List<Integer> list = new ArrayList<Integer>(sampleToType.keySet());
        Collections.sort(list);
        fu.setOutput(targetFileName);
        fu.printLine("SampleID\tLuminal");
        for (Integer sample : list) {
            // Need a sample just before sample id for downstream data processing
            fu.printLine(("Sample " + sample) + "\t" + sampleToType.get(sample));
        }
        fu.close();
    }
    
}
