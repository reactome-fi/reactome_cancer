/*
 * Created on Oct 25, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.r3.fi.SurvivalAnalysisResult;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * Most of methods in this class are used to handle MCL gene expression clusters.
 * @author wgm
 *
 */
public class OvarianMCLGeneExpModuleAnalyzer extends TCGAOvarianCancerAnalyzer {
    public static final int MCL_SIZE_CUTOFF = 8;
    public static final double MCL_MEAN_CUTOFF = 0.25d;
    
    public OvarianMCLGeneExpModuleAnalyzer() {
    }
    
    @Test
    public void doOneStepTrainingAndValidation() throws Exception {
        MCLClusterWrapper mclCluster = new MCLClusterWrapper();
        mclCluster.setUseAbsolute(true);
        
        // Validation with two data sets
        String[] geneExpFileNames = new String[] {
                OVARIAN_DIR_NAME + "TothillDataset/GSE9891_Gene_Exp_111811.txt",
                OVARIAN_DIR_NAME + "data_072610/TCGA_batch9-15_17-19_21-22_24.UE.No_Normal.txt",
                OVARIAN_DIR_NAME + "GSE13876/GSE13876_Exp_log2_norm_gene_avg.txt",
                OVARIAN_DIR_NAME + "GSE26712/GSE26712_Gene_Exp_111911.txt"
        };
        String[] clinFileNames = new String[] {
                OVARIAN_DIR_NAME + "TothillClinInfo111811.txt",
                OVARIAN_DIR_NAME + "data_031910/Batches9-22_tcga_OV_clinical_csv.2010-01-25-noDates.txt",
                OVARIAN_DIR_NAME + "GSE13876/GSE13876_Clin.txt",
                OVARIAN_DIR_NAME + "GSE26712/GSE26712ClinInfo.txt"
        };
        String[] datasetNames = new String[] {
                "GSE9891",
                "TCGA",
                "GSE13876",
                "GSE26712"
        };
        
        int trainingFileIndex = 1;
        
        List<Set<String>> clusters = mclCluster.mclCluster(geneExpFileNames[trainingFileIndex],
                                                           MCL_SIZE_CUTOFF,
                                                           MCL_MEAN_CUTOFF);
//        String clusterFileName = OVARIAN_DIR_NAME + datasetNames[trainingFileIndex] + "MCLClusters_040512.txt";
        String clusterFileName = OVARIAN_DIR_NAME + datasetNames[trainingFileIndex] + "MCLSingleWeightClusters_040512.txt";
        new NetworkClusterAnalyzer().outputNetworkClusters(clusters, clusterFileName);
//        String fileNamePostFix = "SampleToMCLModules_040512.txt";
        String fileNamePostFix = "SampleToSingleWeightModules_040512.txt";
        String tmpFileName = OVARIAN_DIR_NAME + datasetNames[trainingFileIndex] + fileNamePostFix;
        String[] survivalResults = doSurvivalAnalysis(geneExpFileNames[trainingFileIndex], 
                                                      clinFileNames[trainingFileIndex],
                                                      tmpFileName,
                                                      clusters,
                                                      "coxph",
                                                      null);
        System.out.println("Training file: " + geneExpFileNames[trainingFileIndex]);
        System.out.println("Total clusters: " + clusters.size());
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println(i + "\t" + clusters.get(i).size());
        }
        System.out.println("\nSurvival Results from coxph:");
        for (String text : survivalResults)
            System.out.println(text);
        
        for (int i = 0; i < geneExpFileNames.length; i++) {
            if (i == trainingFileIndex)
                continue;
            String validation = geneExpFileNames[i];
            System.out.println("\n" + validation);
            tmpFileName = OVARIAN_DIR_NAME + datasetNames[i] + fileNamePostFix;
            survivalResults = doSurvivalAnalysis(validation, 
                                                 clinFileNames[i],
                                                 tmpFileName,
                                                 clusters,
                                                 "coxph",
                                                 null);
            for (String text : survivalResults)
                System.out.println(text);
        }
    }
    
    private String[] doSurvivalAnalysis(String geneExpFileName,
                                        String clinFileName,
                                        String tmpFileName,
                                        List<Set<String>> clusters,
                                        String model,
                                        Integer moduleIndex) throws Exception {
        SurvivalAnalysisResult result = CancerAnalysisUtilitites.doGeneExpClusterSurvivalAnalysis(geneExpFileName, 
                                                                                                  clinFileName,
                                                                                                  model,
                                                                                                  moduleIndex + "",
                                                                                                  null,
                                                                                                  clusters);
        return new String[]{result.getOutput(), result.getError()};
    }
    
    /**
     * This method is used to do overlapping analysis for selected genes from 
     * R superpc analysis results using script SuperpcClassifier.R.
     * @throws Exception
     */
    @Test
    public void overlapAnalysisForSuperpcGenes() throws Exception {
        String fileName = OVARIAN_DIR_NAME + "SelectedMCLModules.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split(", ");
        Set<String> selectedMCLClusterGenes = new HashSet<String>();
        for (String token : tokens)
            selectedMCLClusterGenes.add(token);
        System.out.println("Total genes in selected MCL clusters: " + selectedMCLClusterGenes.size());
        fu.close();
        fileName = OVARIAN_DIR_NAME + "TCGAOVSuperPCSelectedGenes_1.6.txt";
        Set<String> selectedSingleGenes = fu.loadInteractions(fileName);
        System.out.println("Total genes on gene expression: " + selectedSingleGenes.size());
        // Check overlapping
        // Get total genes
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        selectedSingleGenes.retainAll(totalGenes);
        System.out.println(" in the FI network: " + selectedSingleGenes.size());
        Set<String> shared = InteractionUtilities.getShared(selectedMCLClusterGenes, 
                                                            selectedSingleGenes);
        System.out.println("Shared genes: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(),
                                                                    selectedMCLClusterGenes.size(),
                                                                    selectedSingleGenes.size(), 
                                                                    shared.size());
        System.out.println("P-value: " + pvalue);
        System.out.println("Shared genes are: " + InteractionUtilities.joinStringElements(", ", shared));
    }
    
    /**
     * This method is used to permutate expression correlation for FIs.
     * @throws Exception
     */
    @Test
    public void permutateExpCorrelationAmongFIs() throws Exception {
        String srcFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        List<Double> values = new ArrayList<Double>();
        List<String> fis = new ArrayList<String>();
        fu.setInput(srcFileName);
        int index = 0;
        String line = null;
        while ((line = fu.readLine()) != null) {
            index = line.lastIndexOf("\t");
            values.add(new Double(line.substring(index + 1)));
            fis.add(line.substring(0, index));
        }
        fu.close();
        RandomData randomizer = new RandomDataImpl();
        int[] randomIndices = randomizer.nextPermutation(values.size(), values.size());
        String outFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510_random_102710.txt";
        fu.setOutput(outFileName);
        for (int i = 0; i < fis.size(); i++) {
            index = randomIndices[i];
            //System.out.println(i + " -> " + index);
            fu.printLine(fis.get(i) + "\t" + values.get(index));
        }
        fu.close();
    }
    
    /**
     * This method is used to generate FIs and Cytoscape attribute files based on
     * selected mcl clusters.
     * @throws IOException
     */
    @Test
    public void generateFIsForSelectedMCLClusters() throws IOException {
        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        List<Set<String>> clusters = new CancerGeneExpressionCommon().selectGeneExpClusters(clusterFileName,
                                                                                            fiToCorFileName,
                                                                                            MCL_SIZE_CUTOFF,
                                                                                            MCL_MEAN_CUTOFF);
        System.out.println("Total selected clusters: " + clusters.size());
        Set<String> allClusterGenes = new HashSet<String>();
        for (Set<String> cluster : clusters)
            allClusterGenes.addAll(cluster);
        System.out.println("Total genes in selected clusters: " + allClusterGenes.size());
        // Get all FIs among cluster genes
        Set<String> allFIs = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fis = InteractionUtilities.getFIs(allClusterGenes, allFIs);
        System.out.println("Total FIs among total genes: " + fis.size());
        // Output fis
        String fileName = OVARIAN_DIR_NAME + "FIsInSelectedMCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50_110910.txt";
        fu.saveInteractions(fis, fileName);
        // Generate an cluster attribute
        fileName = OVARIAN_DIR_NAME + "GenesInSelectedMCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50_110910.txt";
        fu.setOutput(fileName);
        fu.printLine("Module (class=java.lang.Integer)");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            for (String gene : cluster) {
                fu.printLine(gene + " = "  + i);
            }
        }
        fu.close();
    }
    
    /**
     * This method is used to generate a matrix for sample to MCL modules exported from FI plug-in.
     * @throws Exception
     */
    @Test
    public void generateSampleToMCLModuleExpressMatrixFromPlugInModules() throws Exception {
//        String plugInFile = OVARIAN_DIR_NAME + "TCGAOV27MCLModules_111811.txt";
////        String outFileName = OVARIAN_DIR_NAME + "GSE9899SampleToTCGAModuleExpression_111811.txt";
////        String outFileName = OVARIAN_DIR_NAME + "GSE9899_No_z_SampleToTCGAModuleExpression_111811.txt";
//        String outFileName = OVARIAN_DIR_NAME + "GSE9891SampleToTCGAModuleExpression_111811.txt";
//        String outFileName = OVARIAN_DIR_NAME + "GSE9891SampleToTCGAModuleExpression_111811.txt";
//        TCGAOvarianCancerClinicalAnalyzer loader = new TCGAOvarianCancerClinicalAnalyzer();
//        Map<String, Map<String, Double>> geneToSampleToValue = loader.loadTothillGeneExpData();
        
        String plugInFile = OVARIAN_DIR_NAME + "GSE13876/GSE13876MCLModules_111811.txt";
        String outFileName = OVARIAN_DIR_NAME + "TCGASampleToGSE13876MCLModuleExp_111811.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = loadGeneExp();
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        List<Set<String>> modules = helper.loadMCLModules(plugInFile);
        helper.generateSampleToGeneExpClusters(geneToSampleToValue, outFileName, modules);
    }
    
    /**
     * This method is used to generate sample to average gene expression value in a cluster. This
     * file can be loaded into R for further analysis.
     * @throws Exception
     */
    @Test
    public void generateSampleToGeneExpClusters() throws Exception {
        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_100510_random_102710_I50.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsInGene_041709_I50_110110.txt";
        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510_random_102710.txt";
        //String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCLClusters102610.txt";
        //String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCLClusters_random_102710_I50.txt";
//        String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCL_random_Clusters_FIsWithGeneExpAbsCorr_102510_I50_110110.txt";
//        String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
//        String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
        String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_I50_112111.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = loadGeneExp();
        helper.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName,
                                               geneToSampleToValue,
                                               MCL_SIZE_CUTOFF,
                                               MCL_MEAN_CUTOFF,
                                               output);
//        List<Set<String>> clusters = helper.selectGeneExpClusters(clusterFileName,
//                                                                  fiToCorFileName,
//                                                                  MCL_SIZE_CUTOFF,
//                                                                  MCL_MEAN_CUTOFF);
//        List<Set<String>> clusters = selectRandomGeneExpClusters(clusterFileName, 27);
//        System.out.println("Total selected clusters: " + clusters.size());
//        for (int i = 0; i < clusters.size(); i++) {
//            Set<String> cluster = clusters.get(i);
//            System.out.println(i + ": " + cluster);
//        }
////        if (true)
////            return;
//        // Get a list of samples
//        List<String> samples = new ArrayList<String>();
//        Map<String, Double> sampleToValue = geneToSampleToValue.values().iterator().next();
//        for (String sample : sampleToValue.keySet()) {
//            if (sample.contains(".11A."))
//                continue; // Escape normal samples
//            samples.add(sample);
//        }
//        fu.setOutput(output);
//        StringBuilder builder = new StringBuilder();
//        builder.append("Sample");
//        for (int i = 0; i < clusters.size(); i++)
//            builder.append("\tModule").append(i);
//        fu.printLine(builder.toString());
//        double total = 0.0d;
//        for (String sample : samples) {
//            builder.setLength(0);
//            String sampleLabel = sample.substring(0, 12).replaceAll("\\.", "-");
//            builder.append(sampleLabel);
//            for (Set<String> cluster : clusters) {
//                total = 0.0d;
//                for (String gene : cluster) {
//                    sampleToValue = geneToSampleToValue.get(gene);
//                    if (sampleToValue != null)
//                        total += sampleToValue.get(sample);
//                }
//                builder.append("\t").append(total / cluster.size());
//            }
//            fu.printLine(builder.toString());
//        }
//        fu.close();
    }
    
    @Test
    public void generateGSE13876SampleToGeneExpClusters() throws IOException {
        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
        //        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_100510_random_102710_I50.txt";
        //        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsInGene_041709_I50_110110.txt";
        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        //        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510_random_102710.txt";
        //String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCLClusters102610.txt";
        //String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCLClusters_random_102710_I50.txt";
        //        String output = OVARIAN_DIR_NAME + "SampleToGeneExpMCL_random_Clusters_FIsWithGeneExpAbsCorr_102510_I50_110110.txt";
//        String output = OVARIAN_DIR_NAME + "GSE13876_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_random_102710_I50_110810.txt";
        String output = OVARIAN_DIR_NAME + "GSE13876_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_z_102410_I50_112111.txt";
        
        
        //        List<Set<String>> clusters = selectRandomGeneExpClusters(clusterFileName, 27);
//        String geneExpFileName = OVARIAN_DIR_NAME + "GSE13876/GSE13876_Exp_log2_norm_gene_avg.txt";
        String geneExpFileName = OVARIAN_DIR_NAME + "GSE13876/GSE13876_Exp_log2_norm_gene_avg_z.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(geneExpFileName);
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        helper.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName, 
                                               geneToSampleToValue,
                                               MCL_SIZE_CUTOFF, 
                                               MCL_MEAN_CUTOFF, 
                                               output);
    }
    
    @Test
    public void generateGSE26712SampleToGeneExpClusters() throws IOException {
        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        
        String geneExpFileName = OVARIAN_DIR_NAME + "GSE26712/GSE26712_z_Gene_Exp_111911.txt";
        String output = OVARIAN_DIR_NAME + "GSE26712_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_z_102510_I50_112111";
        
//        String geneExpFileName = OVARIAN_DIR_NAME + "GSE26712/GSE26712_Gene_Exp_111911.txt";
//        String output = OVARIAN_DIR_NAME + "GSE26712_no_z_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_z_102510_I50_112111.txt";
        
//        String geneExpFileName = OVARIAN_DIR_NAME + "GSE26712/GSE26712_no_normal_z_Gene_Exp_112111.txt";
//        String output = OVARIAN_DIR_NAME + "GSE26712_no_normal_z_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_z_102510_I50_112111.txt";
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(geneExpFileName);
        helper.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName,
                                               geneToSampleToValue,
                                               MCL_SIZE_CUTOFF,
                                               MCL_MEAN_CUTOFF,
                                               output);
    }
    
    @Test
    public void selectGenesFromMCLModules() throws IOException {
        // These indices are copied from R, and should substract 1 in use in Java
//        int[] moduleIndexes = new int[] {18, 9, 24, 13, 10, 26, 22}; // Results based on clinical information: /data_031910/Batches9-22_tcga_OV_clinical_csv.2010-01-25-noDates.txt
        Set<String> genes = getSelectedGenesFromSuperpc();
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        System.out.println("Total genes: " + geneList.size());
        for (String gene : geneList)
            System.out.println(gene);
    }
    
    /**
     * Get genes selected in superpc analysis.
     * @return
     * @throws IOException
     */
    public Set<String> getSelectedGenesFromSuperpc() throws IOException {
//        int[] moduleIndexes = new int[] {18, 9, 26, 8, 13, 20, 19}; // Results based on clinical information: /data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        
//        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
//        List<Set<String>> clusters = helper.selectGeneExpClusters(clusterFileName,
//                                                                  fiToCorFileName, 
//                                                                  MCL_SIZE_CUTOFF,
//                                                                  MCL_MEAN_CUTOFF);
//        int[] moduleIndexes = new int[] {17, 8, 20, 6, 13, 18}; // Results as of April 5, 2012
        int[] moduleIndexes = new int[] {2, 5, 6, 8, 10, 14, 17, 20, 22, 25}; // Got from stepAIC on Oct 15, 2012
        String clusterFileName = OVARIAN_DIR_NAME + "TCGAMCLClusters_040512.txt";
        List<Set<String>> clusters = new NetworkClusterAnalyzer().loadNetworkClusters(clusterFileName);
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println(i + "\t" + clusters.get(i).size());
        }
        Set<String> genes = new HashSet<String>();
//        System.out.println("ID\tModule");
        for (int i : moduleIndexes) {
            Set<String> cluster = clusters.get(i - 1);
            genes.addAll(cluster);
//            for (String gene : cluster)
//                System.out.println(gene + "\t" + i);
        }
        return genes;
    }
    
    /**
     * This method is used to generate sample to gene expression values in clusters for the Tothill dataset.
     * @throws IOException
     */
    @Test
    public void generateTothillSampleToGeneExpClusters() throws IOException {
//        TCGAOvarianCancerClinicalAnalyzer helper = new TCGAOvarianCancerClinicalAnalyzer();
//        Map<String, Map<String, Double>> geneToSampleToExp = helper.loadTothillGeneExpData();
//        Map<String, String> sampleToStatus = helper.loadTothillSampleToStatus();
//        Map<String, Double> sampleToSurvival = helper.loadTothillSampleToDeath();
//        String output = OVARIAN_DIR_NAME + "TothillSampleToGeneExpMCLClusters102610.txt";
        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithTothillGeneExpAbsCorr_I50_110211.txt";
        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithTothillGeneExpAbsCorr_110211.txt";
//        String output = OVARIAN_DIR_NAME + "TothillSampleToGeneExpMCLClusters_random_102710.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_100510_random_102710_I50.txt";
//        String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510_random_102710.txt";
//        String output = OVARIAN_DIR_NAME + "TothillSampleToGeneExpMCL_random_Clusters_FIsWithGeneExpAbsCorr_102510_I50_110110.txt";
//        String output = OVARIAN_DIR_NAME + "TothillSampleToGeneExpMCL_Clusters_FIsWithTothillGeneExpAbsCorr_I50_110211.txt";
        
        String geneExpFileName = OVARIAN_DIR_NAME + "TothillDataset/GSE9891_z_Gene_Exp_111811.txt";
        String output = OVARIAN_DIR_NAME + "TothillSampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50_112111.txt";
        
//        String geneExpFileName = OVARIAN_DIR_NAME + "TothillDataset/GSE9891_Gene_Exp_111811.txt";
//        String output = OVARIAN_DIR_NAME + "Tothill_no_z_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50_112111.txt";
        
//        String geneExpFileName = OVARIAN_DIR_NAME + "TothillDataset/GSE9891_Gene_Exp_111811.txt";
//        String output = OVARIAN_DIR_NAME + "Tothill_no_z_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50_112111.txt";
//        String clusterFileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsInGene_041709_I50_110110.txt";
        
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        
        Map<String, Map<String, Double>> geneToSampleToExp = common.loadGeneExp(geneExpFileName);
        common.generateSampleToGeneExpClusters(clusterFileName, 
                                               fiToCorFileName, 
                                               geneToSampleToExp, 
                                               MCL_SIZE_CUTOFF, 
                                               MCL_MEAN_CUTOFF, 
                                               output);
        
        // This is used to generate clin information file
//        output = OVARIAN_DIR_NAME + "TothillClinInfo102610.txt";
//        builder.setLength(0);
//        builder.append("Sample\tOSDURATION_MONTHS\tOSEVENT");
//        fu.setOutput(output);
//        fu.printLine(builder.toString());
//        for (String sample : samples) {
//            String status = sampleToStatus.get(sample);
//            if (status.equals("D"))
//                status = "1";
//            else
//                status = "0";
//            Double survival = sampleToSurvival.get(sample);
//            builder.setLength(0);
//            builder.append(sample).append("\t").append(survival).append("\t").append(status);
//            fu.printLine(builder.toString());
//        }
//        fu.close();
    }
    
    private List<Set<String>> selectRandomGeneExpClusters(String clusterFileName,
                                                          int clusterNumber) throws IOException {
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
        // Do a size filtering
        for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
            Set<String> set = it.next();
            if (set.size() <= MCL_SIZE_CUTOFF)
                it.remove();
        }
        RandomData randomizer = new RandomDataImpl();
        int[] randomIndices = randomizer.nextPermutation(clusters.size(), clusterNumber);
        List<Set<String>> rtn = new ArrayList<Set<String>>();
        for (int index : randomIndices) {
            System.out.println(index);
            rtn.add(clusters.get(index));
        }
        return rtn;
    }
    
    @Test
    public void checkMCLClusters() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME); 
//        Map<String, Double> fiToValue = loadFIToValue(OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt");
        Map<String, Double> fiToValue = loadFIToValue(OVARIAN_DIR_NAME + "FIsWithTothillGeneExpAbsCorr_110211.txt");
//        String fileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
//        String fileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_100610.txt";
//        String fileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpCorr_I50_110210.txt";
        String fileName = OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithTothillGeneExpAbsCorr_I50_110211.txt";
//      String fiToCorFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(fileName);
        int index = 0;
        // used to calculate some statistical values
        DescriptiveStatistics stat = new DescriptiveStatistics();
        System.out.println("Index\tSize\tMean\tMedian\tMax\tMin\tSD\tGenes");
        for (Set<String> cluster : clusters) {
            if (cluster.size() > 7) {
                Set<String> fisInCluster = InteractionUtilities.getFIs(cluster, fis);
                for (String fi : fisInCluster) {
                    Double value = fiToValue.get(fi);
                    if (value != null)
                        stat.addValue(value);
                }
                System.out.println(index + "\t" +
                                   cluster.size() + "\t" +
                                   stat.getMean() + "\t" +
                                   stat.getPercentile(0.50d) + "\t" +
                                   stat.getMax() + "\t" + 
                                   stat.getMin() + "\t" +
                                   stat.getStandardDeviation() + "\t" +
                                   InteractionUtilities.joinStringElements(", ", cluster));
                stat.clear();
                index ++;
            }
        }
    }
    
    /**
     * This method is used to increase the original expression correlation by the min value so that
     * all weights used by MCL are non-negative.
     * @throws IOException
     */
    @Test
    public void generatePositiveFIsWithGeneExpCorr() throws IOException {
        String srcFileName = OVARIAN_DIR_NAME + "FIsWithGeneExpCorr_100510.txt";
        // Get the min value
        fu.setInput(srcFileName);
        String line = null;
        int index = 0;
        double min = 0.0d;
        while ((line = fu.readLine()) != null) {
            index = line.lastIndexOf("\t");
            Double value = new Double(line.substring(index + 1));
            if (value < min)
                min = value;
        }
        fu.close();
        String fileName = OVARIAN_DIR_NAME + "FIsWithPosGeneExpCorr_110210.txt";
        fu.setInput(srcFileName);
        fu.setOutput(fileName);
        while ((line = fu.readLine()) != null) {
            index = line.lastIndexOf("\t");
            Double value = new Double(line.substring(index + 1));
            fu.printLine(line.substring(0, index) + "\t" + 
                         (value - min));
        }
        fu.close();
    }
    
    /**
     * This method is used to generate a file containing FIs with gene expression correction as
     * weights for FIs. This file will be feed into MCL for clustering.
     * @throws Exception
     */
    @Test
    public void generateFIsWithGeneExpCorr() throws Exception {
        //Map<String, Map<String, Double>> geneToSampleToValue = loadGeneExp();
        TCGAOvarianCancerClinicalAnalyzer helper = new TCGAOvarianCancerClinicalAnalyzer();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadTothillGeneExpData();
        //String output = OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        String output = OVARIAN_DIR_NAME + "FIsWithTothillGeneExpAbsCorr_110211.txt";
        CancerGeneExpressionCommon common = new CancerGeneExpressionCommon();
        common.generateFIsWithGeneExpCorr(geneToSampleToValue, 
                                          output,
                                          ".11A.");
    }
    
    private Map<String, Double> loadFIToValue(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, Double> fiToValue = new HashMap<String, Double>();
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.lastIndexOf("\t");
            String fi = line.substring(0, index);
            String value = line.substring(index + 1);
            fiToValue.put(fi, new Double(value));
        }
        fu.close();
        return fiToValue;
    }
    
    /**
     * This method is used to check fumnd values in the GSE13876 data set. The fumnd values
     * should be overall survival months as described in paper: Survival-related profile, pathways,
     * and trancription factors in ovarian caner, PLos Medicine.
     * @throws IOException
     */
    @Test
    public void checkFumndInGSE13876() throws IOException {
        String fileName = OVARIAN_DIR_NAME + "GSE13876/fumnd.txt";
        fu.setInput(fileName);
        String patientLine = fu.readLine();
        String statusLine = fu.readLine();
        String fumndLine = fu.readLine();
        fu.close();
        // Get patient information
        List<Integer> patientIds = new ArrayList<Integer>();
        String[] tokens = patientLine.split("\t");
        int index = 0;
        for (String token : tokens) {
            if (token.startsWith("!"))
                continue;
            index = token.indexOf(":");
            String id = token.substring(index + 1, token.length() - 1).trim();
            patientIds.add(new Integer(id));
        }
        // Get status information
        List<Integer> status = new ArrayList<Integer>();
        tokens = statusLine.split("\t");
        for (String token : tokens) {
            if (token.startsWith("!"))
                continue;
            index = token.indexOf(":");
            token = token.substring(index + 1, token.length() - 1).trim();
            status.add(new Integer(token));
        }
        // Get fumnd values
        tokens = fumndLine.split("\t");
        List<Integer> fumndValues = new ArrayList<Integer>();
        for (String token : tokens) {
            if (token.startsWith("!"))
                continue;
            index = token.indexOf(":");
            token = token.substring(index + 1, token.length() - 1).trim();
            Integer value = new Integer(token);
            fumndValues.add(value);
        }
        // Create a map
        Map<Integer, Integer> idToValue = new HashMap<Integer, Integer>();
        Map<Integer, Integer> idToStatus = new HashMap<Integer, Integer>();
        for (int i = 0; i < fumndValues.size(); i++) {
            Integer id = patientIds.get(i);
            Integer value = fumndValues.get(i);
            // Check fumnd values
            if (idToValue.containsKey(id)) {
                Integer oldValue = idToValue.get(id);
                if (!oldValue.equals(value))
                    System.out.println(id + " has different values: " + oldValue + " vs " + value);
            }
            else
                idToValue.put(id, value);
            // Check status
            Integer stat = status.get(i);
            if (idToStatus.containsKey(id)) {
                Integer oldStat = idToStatus.get(id);
                if (!oldStat.equals(stat)) {
                    System.err.println(id + " has different status!");
                }
            }
            else
                idToStatus.put(id, stat);
        }
        
        DescriptiveStatistics stat = new DescriptiveStatistics();
        int count = 0;
        int fiveYear = 0;
        int validCount = 0;
        for (Integer id : idToValue.keySet()) {
            Integer value = idToValue.get(id);
            Integer statValue = idToStatus.get(id);
            System.out.println(count + ": " + value + ", " + statValue);
            stat.addValue(value);
            count ++;
            if (statValue == 1) {
                validCount ++;
                if (value >= 60)
                    fiveYear ++;
            }
            else if (value >= 60) {
                fiveYear ++;
                validCount ++;
            }
        }
        System.out.println("Min: " + stat.getMin());
        System.out.println("Max: " + stat.getMax());
        System.out.println("Average: " + stat.getMean());
        System.out.println("Median: " + stat.getPercentile(50.0d));
        System.out.println("Five Year total: " + fiveYear + "(" + ((double)fiveYear / validCount) + ")");
        if (true)
            return;
        // Generate an clinical information file
        String outFileName = OVARIAN_DIR_NAME + "GSE13876/GSE13876_Clin.txt";
        List<Integer> patientIdList = new ArrayList<Integer>(idToValue.keySet());
        Collections.sort(patientIdList);
        fu.setOutput(outFileName);
        // Generate header
        fu.printLine("PatientId\tOSEVENT\tOSDURATION_MONTHS");
        for (Integer id : patientIdList) {
            Integer value = idToValue.get(id);
            Integer statusValue = idToStatus.get(id);
            fu.printLine("X" + id + "\t" + statusValue + "\t" + value);
        }
        fu.close();
    }
    
    /**
     * This method is used to generate a matrix from GSE13876 matrix downloaded from GEO
     * data set.
     * @throws IOException
     */
    @Test
    public void generateGSE13876ExpressionMatrix() throws IOException {
        String dirName = OVARIAN_DIR_NAME  + "GSE13876/";
        String[] fileNames = new String[] {
                dirName + "GSE13876_series_matrix-1.txt",
                dirName + "GSE13876_series_matrix-2.txt"
        };
        String destFileName = dirName + "GSE13876_Exp.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(destFileName);
        String line = null;
        List<String> list1 = new ArrayList<String>();
        List<String> list2 = new ArrayList<String>();
        String sampleLine1 = null;
        String sampleLine2 = null;
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            while ((line = fu.readLine()) != null) {
                line = line.trim();
                if (line.length() == 0 || line.startsWith("\"ID_REF\""))
                    continue;
                if (line.startsWith("!")) {
                    if (line.contains("assigned unique patient id")) {
                        if (sampleLine1 == null)
                            sampleLine1 = line;
                        else
                            sampleLine2 = line;
                    }
                }
                else {
                    if (fileName.endsWith("-1.txt")) {
                        list1.add(line);
                    }
                    else {
                        list2.add(line);
                    }
                }
            }
            fu.close();
        }
        // Get the samples
        String[] tokens = sampleLine1.split("\t");
        int index = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (String token : tokens) {
            if (token.startsWith("!"))
                continue;
            index = token.indexOf(":");
            token = token.substring(index + 1, token.length() - 1).trim();
            builder.append("\t").append(token);
        }
        tokens = sampleLine2.split("\t");
        for (String token : tokens) {
            if (token.startsWith("!"))
                continue;
            index = token.indexOf(":");
            token = token.substring(index + 1, token.length() - 1).trim();
            builder.append("\t").append(token);
        }
        outFu.printLine(builder.toString());
        for (int i = 0; i < list1.size(); i++) {
            String line1 = list1.get(i);
            String line2 = list2.get(i);
            index = line2.indexOf("\t");
            outFu.printLine(line1 + line2.substring(index));
        }
        outFu.close();
    }
    
    /**
     * This method is used to generate gene to express values based on annotation file.
     * @throws Exception
     */
    @Test
    public void generateGeneToExpValuesForGSE13876() throws Exception {
        String dirName = OVARIAN_DIR_NAME + "GSE13876/";
        String annotFileName = dirName + "GSE13876_Annot.txt";
        fu.setInput(annotFileName);
        String line = null;
        Map<String, String> probeIdToGene = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#") || line.startsWith("!") || line.startsWith("ID"))
                continue;
            String[] tokens = line.split("\t");
            if (tokens[6].trim().length() == 0 || tokens[6].contains(","))
                continue;
            probeIdToGene.put(tokens[0], tokens[6]);
        }
        fu.close();
        Map<String, Set<String>> geneToProbeIds = new HashMap<String, Set<String>>();
        for (String id : probeIdToGene.keySet()) {
            String gene = probeIdToGene.get(id);
            InteractionUtilities.addElementToSet(geneToProbeIds, gene, id);
        }
        System.out.println("Total genes: " + geneToProbeIds.size());
        // Load all expression lines
        String srcFileName = dirName + "GSE13876_Exp_log2_norm.txt";
        fu.setInput(srcFileName);
        line = fu.readLine();
        // Get samples
        String[] tokens = line.split("\t");
        Map<String, Set<Integer>> sampleToIndices = new HashMap<String, Set<Integer>>();
        for (int i = 1; i < tokens.length; i++) {
            String sample = tokens[i];
            int index = sample.indexOf(".");
            if (index > 0)
                sample = sample.substring(0, index);
            InteractionUtilities.addElementToSet(sampleToIndices, sample, i);
        }
        List<String> sampleList = new ArrayList<String>(sampleToIndices.keySet());
        Collections.sort(sampleList);
        // Process express values
        Map<String, String[]> idToValue = new HashMap<String, String[]>();
        while ((line = fu.readLine()) != null) { 
            tokens = line.split("\t");
            // Check the probe id
            String gene = probeIdToGene.get(tokens[0]);
            if (gene == null)
                continue; // Escape this line
            idToValue.put(tokens[0], tokens);
        }
        fu.close();
        // Average based on genes
        Map<String, Double[]> geneToValue = new HashMap<String, Double[]>();
        for (String gene : geneToProbeIds.keySet()) {
            Set<String> ids = geneToProbeIds.get(gene);
            String id0 = ids.iterator().next();
            String[] valueText = idToValue.get(id0);
            Double[] values = new Double[valueText.length - 1];
            for (int i = 1; i < valueText.length; i++) {
                values[i - 1] = new Double(valueText[i]);
            }
            for (String id : ids) {
                if (id.equals(id0))
                    continue;
                valueText = idToValue.get(id);
                for (int i = 1; i < valueText.length; i++) {
                    values[i - 1] += new Double(valueText[i]);
                }
            }
            if (ids.size() > 1) {
                for (int i = 0; i < values.length; i++)
                    values[i] /= ids.size();
            }
            // Average based on samples
            Double[] sampleValues = new Double[sampleList.size()];
            for (int i = 0; i < sampleList.size(); i++) {
                String sample = sampleList.get(i);
                Set<Integer> indices = sampleToIndices.get(sample);
                Integer index0 = indices.iterator().next();
                sampleValues[i] = values[index0 - 1];
                for (Integer index : indices) {
                    if (index.equals(index0))
                        continue;
                    sampleValues[i] += values[index - 1];
                }
                sampleValues[i] /= indices.size();
            }
            geneToValue.put(gene, sampleValues);
        }
        List<String> geneList = new ArrayList<String>(geneToValue.keySet());
        Collections.sort(geneList);
        String destFileName = dirName + "GSE13876_Exp_log2_norm_gene_avg.txt";
        fu.setOutput(destFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (String sample : sampleList)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        for (String gene : geneList) {
            builder.setLength(0);
            builder.append(gene);
            Double[] values = geneToValue.get(gene);
            for (Double value : values)
                builder.append("\t").append(value);
            fu.printLine(builder.toString());
                   
        }
        fu.close();
    }
    
    public String[] getSourceFileNames() {
        String[] fileNames = new String[] {
                OVARIAN_DIR_NAME + "data_072610/TCGA_batch9-15_17-19_21-22_24.UE.No_Normal.txt",
                OVARIAN_DIR_NAME + "GSE13876/GSE13876_Exp_log2_norm_gene_avg_z.txt",
                OVARIAN_DIR_NAME + "GSE26712/GSE26712_z_Gene_Exp_111911.txt",
                OVARIAN_DIR_NAME + "TothillDataset/GSE9891_z_Gene_Exp_111811.txt"
        };
        return fileNames;
    }
    
    public String[] getClinFileNames() {
        String[] fileNames = new String[] {
                OVARIAN_DIR_NAME + "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt",
                OVARIAN_DIR_NAME + "GSE13876/GSE13876_Clin.txt",
                OVARIAN_DIR_NAME + "GSE26712/GSE26712ClinInfo.txt",
                OVARIAN_DIR_NAME + "TothillClinInfo111811.txt"
        };
        return fileNames;
    }
    
    public String[] getDataSetNames() {
        String[] datasetNames = new String[] {
                "TCGA",
                "GSE13876",
                "GSE26712",
                "GSE9891"
        };
        return datasetNames;
    }
}
