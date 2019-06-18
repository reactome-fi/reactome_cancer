/*
 * Created on May 29, 2013
 *
 */
package org.reactome.cancer;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.fi.SurvivalAnalysisResult;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to check some results from the FireHose TCGA data sets.
 * @author gwu
 *
 */
public class TCGAFireHoseOVAnalyzer extends TCGAOvarianCancerAnalyzer {
    private final String FH_DIR = R3Constants.OVARIAN_DIR_NAME + "FireHose/052013/";
    
    public TCGAFireHoseOVAnalyzer() {
    }
    
    @Test
    public void doMCLClusteringBasedOnPValue() throws Exception {
        String scoreFileName = FH_DIR + "CNVGeneVsSurvivalWithFDR.txt";
        Map<String, Double> geneToPvalue = loadGeneToValue(scoreFileName, 0, 2);
        String dataFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "all_data_by_genes_transformed.txt";
        String clinFile = FH_DIR + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2013042100.0.0/" + 
                "OV.clin.transformed.txt";
        
        SurvivalAnalysisResult results = doMCLClusteringBasedOnPValue(geneToPvalue,
                                                                      dataFile,
                                                                      clinFile,
                                                                      null,
                                                                      null);
        System.out.println("Error: \n" + results.getError());
        System.out.println("Result: \n" + results.getOutput());
    }

    protected SurvivalAnalysisResult doMCLClusteringBasedOnPValue(Map<String, Double> geneToPvalue,
                                                                  String dataFile,
                                                                  String clinFile,
                                                                  List<String> trainingSamples,
                                                                  List<Set<String>> rtnClusters) throws Exception {
        MCLClusterWrapper mclWrapper = new MCLClusterWrapper();
        Double valueCutoff = 1.3d; // two genes should have p-value <= 0.05
        //        double valueCutoff = 1.0d;
        //        valueCutoff = null;
        int sizeCutoff = 3;
        
        List<Set<String>> clusters = mclWrapper.mclClusterForGeneScores(geneToPvalue,
                                                                        true, 
                                                                        sizeCutoff,
                                                                        valueCutoff);
//        // Do a filtering based on average gene cluster scores
//        for (Iterator<Set<String>> it = clusters.iterator(); it.hasNext();) {
//            Set<String> cluster = it.next();
//            double total = 0.0d;
//            for (String gene : cluster)
//                total += -Math.log10(geneToPvalue.get(gene));
//            if (total / cluster.size() < valueCutoff)
//                it.remove();
//        }
        
        System.out.println("Total clusters: " + clusters.size());
        System.out.println("\nModule\tSize\tAverage_Gene_Score\tGenes");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            double total = 0.0d;
            for (String gene : cluster)
                total += -Math.log10(geneToPvalue.get(gene));
            System.out.println(i + "\t" + cluster.size() + "\t" + total / cluster.size() + "\t" + cluster);
        }
        if (rtnClusters != null) {
            rtnClusters.clear();
            rtnClusters.addAll(clusters);
        }
        SurvivalAnalysisResult results = CancerAnalysisUtilitites.doGeneExpClusterSurvivalAnalysis(dataFile,
                                                                                                   clinFile, 
                                                                                                   "coxph", 
                                                                                                   null, 
                                                                                                   trainingSamples,
                                                                                                   clusters);
        return results;
    }
    
    /**
     * Get a list of genes from three data types: CNVs, methylation, and mutation.
     * @throws IOException
     */
    @Test
    public void generateGeneListFromMultipleSources() throws Exception {
        String cnvCorrFileName = FH_DIR + "CorrelationBetweenCNVsAndmRNAExp.txt";
        Map<String, Double> cnvGeneToCorr = loadGeneToValue(cnvCorrFileName, 1, 2);
        double threshold = 0.60d;
        Set<String> cnvGenes = filterGenesOnCorr(cnvGeneToCorr, threshold, true);
        System.out.println("CNV genes with corr cutoff " + threshold + ": "  + cnvGenes.size());
        
        String metCorrFileName = FH_DIR + "CorrelationBetweenMaxMethylationAndmRNAExp.txt";
        Map<String, Double> metGeneToCorr = loadGeneToValue(metCorrFileName, 1, 2);
        Set<String> methGenes = filterGenesOnCorr(metGeneToCorr, -0.50d, false);
        System.out.println("Methylation genes with corr cutoff -0.50d: " + methGenes.size());
        
        MATFileLoader matFileLoader = new MATFileLoader();
        String mutationFileName = FH_DIR + "gdac.broadinstitute.org_OV-TP.Mutation_Assessor.Level_4.2013032600.0.0/OV-TP.maf.annotated";
        Map<String, Set<String>> sampleToGenes = matFileLoader.loadSampleToGenes(mutationFileName, false);
        // Check genes mutated at least in 4 or more samples
        int sampleCutoff = 4;
        filterSampleToGenes(sampleToGenes, sampleCutoff);
        Set<String> mutatedGenes = InteractionUtilities.grepAllGenes(sampleToGenes);
        System.out.println("Mutated genes with sample cutoff " + sampleCutoff + ": " + mutatedGenes.size());
        
        // Check gene overlapping
        Set<String> shared = InteractionUtilities.getShared(cnvGenes, mutatedGenes);
        System.out.println("Shared between CNVs and mutated: " + shared.size());
        double pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                                    mutatedGenes.size(),
                                                                    cnvGenes.size(),
                                                                    shared.size());
        System.out.println("p-value: " + pvalue);
        
        shared = InteractionUtilities.getShared(cnvGenes, methGenes);
        System.out.println("Shared between CNVs and Methylation: " + shared.size());
        pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                             cnvGenes.size(),
                                                             methGenes.size(),
                                                             shared.size());
        System.out.println("P-value: " + pvalue);
        
        shared = InteractionUtilities.getShared(mutatedGenes, methGenes);
        System.out.println("Shared between Mutated and Methylation: " + shared.size());
        pvalue = MathUtilities.calculateHypergeometricPValue(R3Constants.TOTAL_HUMAN_GENES,
                                                             mutatedGenes.size(),
                                                             methGenes.size(),
                                                             shared.size());
        System.out.println("P-value: " + pvalue);
    }
    
    /**
     * Filter out genes based a correlation threshold.
     * @param geneToCorr
     * @param corrTheshold
     * @return
     */
    private Set<String> filterGenesOnCorr(Map<String, Double> geneToCorr,
                                          Double corrTheshold,
                                          boolean greater)  {
        Set<String> rtn = new HashSet<String>();
        for (String gene : geneToCorr.keySet()) {
            Double corr = geneToCorr.get(gene);
            if (greater && corr >= corrTheshold)
                rtn.add(gene);
            else if (!greater && corr <= corrTheshold)
                rtn.add(gene);
        }
        return rtn;
    }
    
    /**
     * Load a mapping file from genes to their correlations calculated in some data source.
     * @param fileName
     * @return
     * @throws IOException
     */
    private Map<String, Double> loadGeneToValue(String fileName,
                                                int geneColIndex,
                                                int scoreColIndex) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Double> geneToCorr = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[geneColIndex];
            Double corr = new Double(tokens[scoreColIndex]);
            geneToCorr.put(gene, corr);
        }
        fu.close();
        return geneToCorr;
    }
    
    @Test
    public void checkMethyatlionAndRNACorrelation() throws Exception {
        CancerGeneExpressionCommon loader = new CancerGeneExpressionCommon();
        String geneExpFile = FH_DIR + "gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_4.2013042100.0.0/" + 
                             "OV.uncv2.mRNAseq_RSEM_normalized_log2.transformed.txt";
        Map<String, Map<String, Double>> expGeneToSampleToValue = loader.loadGeneExp(geneExpFile, false);
        System.out.println("Total genes in the gene expression file: " + expGeneToSampleToValue.size());
        if (expGeneToSampleToValue.keySet().contains("NR5A1"))
            System.out.println("Correct in gene exp");
        String methylationFile = FH_DIR + "methylation_jhu_usc_.Level_3.2013042100.0.0/" +
                                 "OV.methylation.max.transformed.txt";
        Map<String, Map<String, Double>> methylationGeneToSampleToValue = loader.loadGeneExp(methylationFile, false);
        System.out.println("Total genes in the methylation file: " + methylationGeneToSampleToValue.size());
        if (methylationGeneToSampleToValue.keySet().contains("NR5A1"))
            System.out.println("Correct in methyation!");
        // Get shared genes in these two files
        Set<String> sharedGenes = new HashSet<String>(expGeneToSampleToValue.keySet());
        sharedGenes.retainAll(methylationGeneToSampleToValue.keySet());
        System.out.println("Total shared genes: " + sharedGenes.size());
        if (sharedGenes.contains("NR5A1"))
            System.out.println("Shared genes has NR5A1");
        final Map<String, Double> geneToCorr = new HashMap<String, Double>();
        List<Double> list1 = new ArrayList<Double>();
        List<Double> list2 = new ArrayList<Double>();
        for (String gene : sharedGenes) {
            Map<String, Double> expSampleToValue = expGeneToSampleToValue.get(gene);
            Map<String, Double> cnvSampleToValue = methylationGeneToSampleToValue.get(gene);
            for (String sample : expSampleToValue.keySet()) {
                Double cnv = cnvSampleToValue.get(sample);
                if (cnv == null)
                    continue;
                list1.add(expSampleToValue.get(sample));
                list2.add(cnv);
            }
            SpearmansCorrelation corr = MathUtilities.constructSpearmansCorrelation(list1, list2);
            geneToCorr.put(gene, corr.getCorrelationMatrix().getEntry(0, 1));
            list1.clear();
            list2.clear();
        }
        List<String> geneList = new ArrayList<String>(geneToCorr.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double corr1 = geneToCorr.get(gene1);
                Double corr2 = geneToCorr.get(gene2);
                return corr1.compareTo(corr2);
            }
        });
        int index = 1;
        String outFile = FH_DIR + "CorrelationBetweenMaxMethylationAndmRNAExp.txt";
//        String outFile = FH_DIR + "CorrelationBetweenMaxMethylationAndmRNAExp.txt";
        fu.setOutput(outFile);
        fu.printLine("Order\tGene\tCorr");
        for (String gene : geneList)
            fu.printLine(index ++ + "\t" + gene + "\t" + geneToCorr.get(gene));
        fu.close();
    }
    
    /**
     * This method is used to merge two gene-based file together.
     * @throws IOException
     */
    @Test
    public void mergeGeneBasedFiles() throws IOException {
//        String[] fileNames = new String[] {
//                FH_DIR + "CNVGeneVsSurvivalWithFDR.txt",
//                FH_DIR + "CNVGeneVariance.txt",
//                FH_DIR + "CorrelationBetweenCNVsAndmRNAExp.txt"
//        };
//        String outFile = FH_DIR + "CNVGenesSurvivalVarianceCorrelation.txt";
        String[] fileNames = new String[] {
                FH_DIR + "ArrayExpGeneVsSurvivalWithFDR.txt",
                FH_DIR + "MethylationGeneVsSurvivalWithFDR.txt",
                FH_DIR + "CNVGeneVsSurvivalWithFDR.txt"
        };
        String outFile = FH_DIR + "GeneBasedSurvivalWithFDR.txt";
        
        List<String> headers = new ArrayList<String>();
        List<Map<String, String>> geneToValueList = new ArrayList<Map<String,String>>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            String line = fu.readLine();
            int index = line.indexOf("\t");
            String header = line.substring(index + 1);
            Map<String, String> geneToValue = new HashMap<String, String>();
            while ((line = fu.readLine()) != null) {
                index = line.indexOf("\t");
                geneToValue.put(line.substring(0, index),
                                line.substring(index + 1));
            }
            fu.close();
            headers.add(header);
            geneToValueList.add(geneToValue);
        }
        
        // Output
        fu.setOutput(outFile);
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        for (String header : headers)
            builder.append("\t").append(header);
        fu.printLine(builder.toString());
        builder.setLength(0);
        boolean shouldExport = true;
        for (String gene : geneToValueList.get(0).keySet()) {
            shouldExport = true;
            builder.setLength(0);
            builder.append(gene);
            for (Map<String, String> geneToValue : geneToValueList) {
                String value = geneToValue.get(gene);
                if (value == null) {
                    shouldExport = false;
                    break;
                }
                builder.append("\t").append(value);
            }
            if (shouldExport)
                fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    /**
     * Use this method to calculate FDRs based on Benjamini Hochberg method.
     * @throws IOException
     */
    @Test
    public void calculateFDRs() throws IOException {
//        String fileName = FH_DIR + "CNVGeneVsSurvival.txt";
//        String outFileName = FH_DIR + "CNVGeneVsSurvivalWithFDR.txt";
        
//        String fileName = FH_DIR + "MethylationGeneVsSurvival.txt";
//        String outFileName = FH_DIR + "MethylationGeneVsSurvivalWithFDR.txt";
        
        String fileName = FH_DIR + "ArrayExpGeneVsSurvival.txt";
        String outFileName = FH_DIR + "ArrayExpGeneVsSurvivalWithFDR.txt";
        
        fu.setInput(fileName);
        String line = fu.readLine();
        List<Double> pvalues = new ArrayList<Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            pvalues.add(new Double(tokens[2]));
        }
        fu.close();
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        // Output
        int index = 0;
        fu.setOutput(outFileName);
        fu.setInput(fileName);
        line = fu.readLine();
        fu.printLine(line + "\tFDR");
        while ((line = fu.readLine()) != null) {
            Double fdr = fdrs.get(index);
            fu.printLine(line + "\t" + fdr);
            index ++;
        }
        fu.close();
    }
    
    @Test
    public void doMCLBasedCrossValidation() throws Exception {
        String cnvFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "focal_data_by_genes_transformed.txt";
        // CNV file should have the same format as gene expression. So we can use the same loader
        // to load this file.
        Map<String, Map<String, Double>> cnvGeneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(cnvFile);
        System.out.println("Total genes in the CNV file: " + cnvGeneToSampleToValue.size());
        Set<String> samples = CancerAnalysisUtilitites.grepSamplesInGeneToSampleToValue(cnvGeneToSampleToValue);
        System.out.println("Total numbers in CNVs: " + samples.size());
        List<Set<String>> cvSamples = CancerAnalysisUtilitites.randomCrossValidationSplit(samples, 3);
        // Check how many samples in the clinical file
        String clinFile = FH_DIR + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2013042100.0.0/" + 
                "OV.clin.transformed.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(cnvFile);
        String outFile = FH_DIR + "CNVSurivlva4FoldCV_MCLCluster_GeneFilter_060413.txt";
        fu.setOutput(outFile);
        for (int i = 0; i < cvSamples.size(); i++) {
            Set<String> validationSamples = cvSamples.get(i);
            fu.printLine("\nValidation fold: " + i + ", " + validationSamples.size() + " samples");
            // Just a simpler way to get the training samples, instead of merging all remaining
            // folds.
            Set<String> trainingSamples = new HashSet<String>(samples);
            trainingSamples.removeAll(validationSamples);
            fu.printLine("Training samples: " + trainingSamples.size());
            fu.printLine("Training:");
            SurvivalAnalysisResult result = checkGenesSurvival(geneToSampleToValue, 
                                                              clinFile, 
                                                              trainingSamples);
            Map<String, Double> geneToPvalue = extractGeneToPValue(result);
            List<Set<String>> rtnClusters = new ArrayList<Set<String>>(); // To hold clustering results
            SurvivalAnalysisResult mclSurvivalResult = doMCLClusteringBasedOnPValue(geneToPvalue, 
                                                                                    cnvFile,
                                                                                    clinFile,
                                                                                    new ArrayList<String>(trainingSamples),
                                                                                    rtnClusters);
            fu.printLine(mclSurvivalResult.getOutput());
            fu.printLine("Validation:");
            result = CancerAnalysisUtilitites.doGeneExpClusterSurvivalAnalysis(cnvFile,
                                                                               clinFile, 
                                                                               "coxph", 
                                                                               null, 
                                                                               new ArrayList<String>(validationSamples),
                                                                               rtnClusters);
            fu.printLine(result.getOutput());
//            break;
        }
        fu.close();
    }
    
    private Map<String, Double> extractGeneToPValue(SurvivalAnalysisResult result) throws IOException {
        Map<String, Double> geneToValue = new HashMap<String, Double>();
        String output = result.getOutput();
        StringReader stringReader = new StringReader(output);
        BufferedReader br = new BufferedReader(stringReader);
        String line = br.readLine(); // Escape the title line
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToValue.put(tokens[0], new Double(tokens[2]));
        }
        br.close();
        return geneToValue;
    }
    
    @Test
    public void doGeneBasedCrossValidation() throws Exception {
        String cnvFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "focal_data_by_genes_transformed.txt";
        // CNV file should have the same format as gene expression. So we can use the same loader
        // to load this file.
        Map<String, Map<String, Double>> cnvGeneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(cnvFile);
        System.out.println("Total genes in the CNV file: " + cnvGeneToSampleToValue.size());
        Set<String> samples = CancerAnalysisUtilitites.grepSamplesInGeneToSampleToValue(cnvGeneToSampleToValue);
        System.out.println("Total numbers in CNVs: " + samples.size());
        List<Set<String>> cvSamples = CancerAnalysisUtilitites.randomCrossValidationSplit(samples, 3);
        // Check how many samples in the clinical file
        String clinFile = FH_DIR + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2013042100.0.0/" + 
                "OV.clin.transformed.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(cnvFile);
        String outFile = FH_DIR + "CNVSurivlva4FoldCV_060413.txt";
        fu.setOutput(outFile);
        for (int i = 0; i < cvSamples.size(); i++) {
            Set<String> validationSamples = cvSamples.get(i);
            fu.printLine("\nValidation fold: " + i + ", " + validationSamples.size() + " samples");
            // Just a simpler way to get the training samples, instead of merging all remaining
            // folds.
            Set<String> trainingSamples = new HashSet<String>(samples);
            trainingSamples.removeAll(validationSamples);
            fu.printLine("Training samples: " + trainingSamples.size());
            fu.printLine("Training:");
            SurvivalAnalysisResult result = checkGenesSurvival(geneToSampleToValue, 
                                                              clinFile, 
                                                              trainingSamples);
            fu.printLine(result.getOutput());
            fu.printLine("Validation:");
            result = checkGenesSurvival(geneToSampleToValue, 
                                        clinFile, 
                                        validationSamples);
            fu.printLine(result.getOutput());
        }
        fu.close();
    }
    
    /**
     * Choose genes that are significant across both training and validation samples from
     * 3-fold gene-based cross-validation.
     * @throws Exception
     */
    @Test
    public void selectCNVSurvivalGenes() throws Exception {
        double pvalueCutoff = 0.10d;
        //String fileName = FH_DIR + "CNVSurivlva3FoldCV_060413_merged.txt";
        String fileName = FH_DIR + "CNVSurivlva4FoldCV_060313_merged.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        int total = 0;
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            String[] tokens = line.split("\t");
            boolean isPicked = true;
            // Make sure all p-values less than the threshold value
            for (int i = 0; i < tokens.length; i++) {
                if (i % 3 == 2) {
                    double pvalue = Double.parseDouble(tokens[i]);
                    if (pvalue > pvalueCutoff) {
                        isPicked = false;
                        break;
                    }
                }
            }
            if (isPicked) {
                System.out.println(tokens[0]);
                total ++;
            }
        }
        fu.close();
        System.out.println("Total picked genes: " + total);
    }
    
    @Test
    public void processCVSurvivalFile() throws Exception {
        String srcFileName = FH_DIR + "CNVSurivlva4FoldCV_060313.txt";
        String outFileName = FH_DIR + "CNVSurivlva4FoldCV_060313_merged.txt";
        
//        String srcFileName = FH_DIR + "CNVSurivlva4FoldCV_MCLCluster_060313.txt";
//        String outFileName = FH_DIR + "CNVSurivlva4FoldCV_MCLCluster_060313_merged.txt";
        
//        String srcFileName = FH_DIR + "CNVSurivlva3FoldCV_060413.txt";
//        String outFileName = FH_DIR + "CNVSurivlva3FoldCV_060413_merged.txt";
        
        fu.setInput(srcFileName);
        fu.setOutput(outFileName);
        String line = null;
        String header = null;
        Map<String, String> geneToLine = new HashMap<String, String>();
        boolean isInData = false;
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0 || line.startsWith("Training") || line.startsWith("Validation")) {
                isInData = false;
                continue;
            }
            if (line.startsWith("Module")) {
                isInData = true;
                if (header == null)
                    header = line;
                else
                    header += line;
                continue;
            }
            if (isInData) {
                String[] tokens = line.split("\t");
                String gene = tokens[0];
                String geneLine = geneToLine.get(gene);
                if (geneLine == null) {
                    geneToLine.put(gene, line);
                }
                else {
                    geneToLine.put(gene, geneLine + line);
                }
            }
        }
        List<String> geneList = new ArrayList<String>(geneToLine.keySet());
        Collections.sort(geneList);
        fu.printLine(header);
        for (String gene : geneList) {
            fu.printLine(geneToLine.get(gene));
        }
        fu.close();
    }
    
    /**
     * This method is used to do a gene-based survival analysis based on CNVs using R.
     * @throws Exception
     */
    @Test
    public void checkGenesSurvival() throws Exception {
        // Gene-based survival analysis for CNVs
//        String dataFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
//                "all_data_by_genes_transformed.txt";
        // Gene-based survival analysis for methylation
//        String dataFile = FH_DIR + "methylation_jhu_usc_.Level_3.2013042100.0.0/" +
//                "OV.methylation.max.transformed.txt";       
        // Gene-based survival analysis for mRNA gene expression
        String dataFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.mRNA_Preprocess_Median.Level_4.2013032600.0.0/" + 
                "OV-TP.medianexp.transformed.txt";
        
        String clinFile = FH_DIR + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2013042100.0.0/" + 
                "OV.clin.transformed.txt";
        
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(dataFile, false);
        SurvivalAnalysisResult result = checkGenesSurvival(geneToSampleToValue, clinFile, null);
        if (result.getError() != null)
            System.out.println("Error: " + result.getError());
        if (result.getOutput() != null)
            System.out.println("Survival analysis results: \n" + result.getOutput());
    }

    private SurvivalAnalysisResult checkGenesSurvival(Map<String, Map<String, Double>> geneToSampleToValue,
                                                      String clinFile,
                                                      Collection<String> samples) throws IOException {
        CancerGeneExpressionCommon loader = new CancerGeneExpressionCommon();
        // Create a temp score file from sample to gene to value
        String scoreFileName = R3Constants.TEMP_DIR + "score.txt";
        if (samples == null)
            loader.generateSampleToGeneToExpValue(geneToSampleToValue, scoreFileName);
        else
            loader.generateSampleToGeneToExpValue(geneToSampleToValue, samples, scoreFileName);
        
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        SurvivalAnalysisResult result = survivalHelper.doSurvivalAnalysis(new File(scoreFileName),
                                                                          new File(clinFile),
                                                                          "coxph",
                                                                          null,
                                                                          "cnv",
                                                                          false);
        return result;
    }
    
    /**
     * This method is used to caculate variations among samples for a data set. The calculation
     * is gene-based.
     * @throws Exception
     */
    @Test
    public void checkVariance() throws Exception {
        String dataFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "all_data_by_genes_transformed.txt";
        String outFileName = FH_DIR + "CNVGeneVariance.txt";
        fu.setOutput(outFileName);
        
        CancerGeneExpressionCommon loader = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = loader.loadGeneExp(dataFile, false);
        SummaryStatistics stat = new SummaryStatistics();
        fu.printLine("Gene\tVariance\tSD");
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            for (String sample : sampleToValue.keySet()) {
                Double value = sampleToValue.get(sample);
                stat.addValue(value);
            }
            fu.printLine(gene + "\t" + stat.getVariance() + "\t" + stat.getStandardDeviation());
            stat.clear();
        }
        fu.close();
    }
    
    /**
     * A simple method to check how many samples in the CNV data file.
     * @throws Exception
     */
    @Test
    public void checkSamplesInCNV() throws Exception {
        String cnvFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "focal_data_by_genes_transformed.txt";
        // CNV file should have the same format as gene expression. So we can use the same loader
        // to load this file.
        Map<String, Map<String, Double>> cnvGeneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(cnvFile);
        System.out.println("Total genes in the CNV file: " + cnvGeneToSampleToValue.size());
        Set<String> samples = CancerAnalysisUtilitites.grepSamplesInGeneToSampleToValue(cnvGeneToSampleToValue);
        System.out.println("Total numbers in CNVs: " + samples.size());
        // Check how many samples in the clinical file
        String clinFile = FH_DIR + "gdac.broadinstitute.org_OV.Clinical_Pick_Tier1.Level_4.2013042100.0.0/" + 
                "OV.clin.transformed.txt";
        Set<String> osSamples = new HashSet<String>();
        fu.setInput(clinFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[0];
            if (sample.endsWith(".1"))
                continue; // Some bugs in the original clinical file
            String osduration = tokens[22];
            if (osduration.equals("NA"))
                continue;
            osSamples.add(sample);
        }
        fu.close();
        System.out.println("Total OS samples: " + osSamples.size());
        samples.retainAll(osSamples);
        System.out.println("CNV samples with OS values: " + samples.size());
    }
    
    /**
     * Check correlations between genes' CNVs and mRNA expressions.
     * @throws Exception
     */
    @Test
    public void checkCNVAndmRNACorrelation() throws Exception {
        CancerGeneExpressionCommon loader = new CancerGeneExpressionCommon();
        String geneExpFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.mRNA_Preprocess_Median.Level_4.2013032600.0.0/" + 
                             "OV-TP.medianexp.transformed.txt";
        Map<String, Map<String, Double>> expGeneToSampleToValue = loader.loadGeneExp(geneExpFile);
        System.out.println("Total genes in the gene expression file: " + expGeneToSampleToValue.size());
        
        // Note: this file, not the threshold one
//        String cnvFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
//                         "all_data_by_genes_transformed.txt";
        String cnvFile = FH_DIR + "gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2013032600.0.0/" +
                "focal_data_by_genes_transformed.txt";
        // CNV file should have the same format as gene expression. So we can use the same loader
        // to load this file.
        Map<String, Map<String, Double>> cnvGeneToSampleToValue = loader.loadGeneExp(cnvFile);
        System.out.println("Total genes in the CNV file: " + cnvGeneToSampleToValue.size());
        
        // Get shared genes in these two files
        Set<String> sharedGenes = new HashSet<String>(expGeneToSampleToValue.keySet());
        sharedGenes.retainAll(cnvGeneToSampleToValue.keySet());
        System.out.println("Total shared genes: " + sharedGenes.size());
        
        // Calculate Pearson correlation
        final Map<String, Double> geneToCorr = new HashMap<String, Double>();
        Map<String, Double> geneToPValue = new HashMap<String, Double>();
        List<Double> list1 = new ArrayList<Double>();
        List<Double> list2 = new ArrayList<Double>();
        for (String gene : sharedGenes) {
            Map<String, Double> expSampleToValue = expGeneToSampleToValue.get(gene);
            Map<String, Double> cnvSampleToValue = cnvGeneToSampleToValue.get(gene);
            for (String sample : expSampleToValue.keySet()) {
                Double cnv = cnvSampleToValue.get(sample);
                if (cnv == null)
                    continue;
                list1.add(expSampleToValue.get(sample));
                list2.add(cnv);
            }
            PearsonsCorrelation corr = MathUtilities.constructPearsonCorrelation(list1, list2);
            geneToCorr.put(gene, corr.getCorrelationMatrix().getEntry(0, 1));
            geneToPValue.put(gene, corr.getCorrelationPValues().getEntry(0, 1));
            list1.clear();
            list2.clear();
        }
        List<String> geneList = new ArrayList<String>(geneToCorr.keySet());
        Collections.sort(geneList, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double corr1 = geneToCorr.get(gene1);
                Double corr2 = geneToCorr.get(gene2);
                return corr2.compareTo(corr1);
            }
        });
        int index = 1;
        String outFile = FH_DIR + "CorrelationBetweenFocalCNVsAndmRNAExp.txt";
        fu.setOutput(outFile);
        fu.printLine("Order\tGene\tCorr\tPValue");
        for (String gene : geneList)
            fu.printLine(index ++ + "\t" + gene + "\t" + geneToCorr.get(gene) + "\t" + geneToPValue.get(gene));
        fu.close();
    }
    
}
