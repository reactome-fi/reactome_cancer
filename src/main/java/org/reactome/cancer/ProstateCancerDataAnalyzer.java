/*
 * Created on Feb 22, 2013
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

import org.apache.commons.math.MathException;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to do some analyses for prostate cancer.
 * @author gwu
 *
 */
public class ProstateCancerDataAnalyzer {
    private final String DIR_NAME = "/Users/gwu/Documents/wgm/work/Keli/";
    private final String TCGA_DIR = "datasets/TCGA/PRAD/gdac.broadinstitute.org_PRAD.mRNAseq_Preprocess.Level_4.2013111400.0.0/";
    private GSESoftTextDataHandler dataHandler = new GSESoftTextDataHandler();
    private FileUtility fu = new FileUtility();
    
    public ProstateCancerDataAnalyzer() {
    }
    
    @Test
    public void checkLFNGExpressionFromGSE() throws Exception {
//        String fileName = DIR_NAME + "GDS2547_LFNG_GSE6919.txt";
        String fileName = DIR_NAME + "GDS2545_LFNG_GSE6919.txt";
        fu.setInput(fileName);
        List<Double> cancer = new ArrayList<Double>();
        List<Double> normal = new ArrayList<Double>();
        List<Double> meta = new ArrayList<Double>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[1].startsWith("Normal prostate tissue adjacent to tumor")) {
                normal.add(new Double(tokens[2]));
            }
            else if (tokens[1].startsWith("Tumor samples")) {
                cancer.add(new Double(tokens[2]));
            }
            else if (tokens[1].startsWith("Metastatic prostate tumor samples")) {
                meta.add(new Double(tokens[2]));
            }
        }
        fu.close();
        System.out.println("Normal: " + normal.size());
        calculateMean(normal);
        System.out.println("Cancer: " + cancer.size());
        calculateMean(cancer);
        System.out.println("Meta: " + meta.size());
        calculateMean(meta);
        System.out.println("T-test between normal and primary:");
        double pvalue = MathUtilities.calculateTTest(normal, cancer);
        System.out.println("P-value: " + pvalue);
        System.out.println("T-test between normal and meta:");
        pvalue = MathUtilities.calculateTTest(normal, meta);
        System.out.println("P-value: " + pvalue);
    }
    
    @Test
    public void outputMSKCCLFNGNKX31AndGleason() throws Exception {
        String fileName = DIR_NAME + "MSKCC/MSKCC_PCa_Clinical_Annotation.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Double> sampleToGleason = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Get primary only
            if (!tokens[1].equals("PRIMARY"))
                continue;
            Double gleason = null;
            try {
                gleason = new Double(tokens[25]);
            }
            catch(NumberFormatException e) {
                e.printStackTrace();
            }
            sampleToGleason.put(tokens[0], gleason); // 25 for PathGGS
        }
        fu.close();
        
//        fileName = DIR_NAME + "MSKCC/MSKCC_PCa_mRNA_Zscores_Modified.txt";
        fileName = DIR_NAME + "MSKCC/MSKCC_PCa_mRNA_data_modified.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName);
        Map<String, Double> lfngExp = geneToSampleToValue.get("LFNG");
        lfngExp = filterMSKCCGeneExp(lfngExp, sampleToGleason.keySet());
        Map<String, Double> nkx31Exp = geneToSampleToValue.get("NKX3-1");
        nkx31Exp = filterMSKCCGeneExp(nkx31Exp, sampleToGleason.keySet());

        System.out.println("Size of gleason: " + sampleToGleason.size());
        System.out.println("Size of LFNG: " + lfngExp.size());
        System.out.println("Size of NKX3-1: " + nkx31Exp.size());
//        if (true)
//            return;
//        String destFileName = DIR_NAME + "MSKCC_LFNG_NKX31_Gleason_122313.txt";
        String destFileName = DIR_NAME + "MSKCC_LFNG_NKX31_Gleason_normal_122913.txt";
        outputLFNGNKX31Gleason(lfngExp, 
                               nkx31Exp,
                               sampleToGleason,
                               destFileName);
    }
    
    private Map<String, Double> filterMSKCCGeneExp(Map<String, Double> sampleToValue,
                                                   Set<String> primarySamples) {
        Map<String, Double> rtn = new HashMap<String, Double>();
        // Want to get primary and normal samples only
        for (String sample : sampleToValue.keySet()) {
             if (sample.startsWith("PAN"))
                 rtn.put(sample + ".normal", sampleToValue.get(sample));
             else if (sample.startsWith("PCA") && primarySamples.contains(sample)) {
                 rtn.put(sample + ".primary", sampleToValue.get(sample));
             }
        }
        // Make sure all primary samples have been entered
        Set<String> copy = new HashSet<String>(primarySamples);
        copy.removeAll(sampleToValue.keySet());
        for (String sample : copy)
            rtn.put(sample + ".primary", null);
        return rtn;
    }
    
    @Test
    public void mergeMSKCCLFNGExpressionAndClinFile() throws Exception {
        String fileName = DIR_NAME + "MSKCC/MSKCC_PCa_mRNA_Zscores_Modified.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName);
        Map<String, Double> lfngSampleToValue = geneToSampleToValue.get("LFNG");
        fileName = DIR_NAME + "MSKCC/MSKCC_PCa_Clinical_Annotation.txt";
        fu.setInput(fileName);
        String destFileName = DIR_NAME + "MSKCC/MSKCC_PCa_Clinical_Annotation_LFNG.txt";
        fu.setOutput(destFileName);
        String line = fu.readLine();
        fu.printLine(line + "\tLFNG");
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                continue;
            String[] tokens = line.split("\t");
            Double value = lfngSampleToValue.get(tokens[0]);
            if (value == null)
                fu.printLine(line + "\tNA");
            else
                fu.printLine(line + "\t" + value);
        }
        fu.close();
    }
    
    @Test
    public void mergeTCGALFNGExpressionAndClinFile() throws Exception {
        String fileName = DIR_NAME + "TCGA_LFNG_Exp.txt";
        Map<String, Double> sampleToExp = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line1 = fu.readLine();
        String[] tokens1 = line1.split("\t");
        String line2 = fu.readLine();
        String[] tokens2 = line2.split("\t");
        for (int i = 1; i < tokens1.length; i++) {
            if (tokens1[i].endsWith("-11"))
                continue; // Normal samples
            String sample = tokens1[i].substring(0, 12).toLowerCase(); // Get rid of the last -01
            Double value = new Double(tokens2[i]);
            sampleToExp.put(sample, value);
        }
        fu.close();
        fileName = "datasets/TCGA/PRAD/gdac.broadinstitute.org_PRAD.Clinical_Pick_Tier1.Level_4.2013111400.0.0/" +
                   "PRAD.clin.merged.picked.transformed.txt";
        String destFileName = DIR_NAME + "TCGA_Clin_LFNG.txt";
        fu.setOutput(destFileName);
        fu.setInput(fileName);
        String line = fu.readLine();
        fu.printLine(line + "\t" + "LFNG");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double value = sampleToExp.get(tokens[0]);
            if (value == null)
                fu.printLine(line + "\tNA");
            else
                fu.printLine(line + "\t" + value);
        }
        fu.close();
    }
    
    @Test
    public void checkGSE6919LFNGGleasonScore() throws Exception {
        List<String> primarySamples = loadGSE6919PrimarySamples();
        Map<String, String> sampleToGleason = loadSampleToFetaure("Gleason Grade");
//        String fileName = DIR_NAME + "GSE6919_GPL8300_121813.txt";
        String fileName = DIR_NAME + "GSE6919_GPL93_121813.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName);
        Map<String, Double> sampleToValue = geneToSampleToValue.get("LFNG");
        List<Double> lfng = new ArrayList<Double>();
        List<Double> gleason = new ArrayList<Double>();
        for (String sample : primarySamples) {
            Double value = sampleToValue.get(sample);
            if (value == null)
                continue;
            String score = sampleToGleason.get(sample);
            if (score == null || score.equals("NA"))
                continue;
            lfng.add(value);
            gleason.add(new Double(score));
        }
        PearsonsCorrelation cor = MathUtilities.constructPearsonCorrelation(lfng, gleason);
        System.out.println("Source file: " + fileName);
        System.out.println("Correlation between LFNG and Gleason Score for " + lfng.size() + " data points:");
        System.out.println(cor.getCorrelationMatrix().getEntry(0, 1) + " with pvalue " + 
                           cor.getCorrelationPValues().getEntry(0, 1));
    }
    
    @Test
    public void checkGSE6919LFNGGeneCorrelations() throws Exception {
        String fileName = DIR_NAME + "GSE6919_Samples_121913.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get a list of primary cancers
        List<String> primary = new ArrayList<String>();
        Map<String, String> sampleToTissue = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[3].equals("primary"))
                primary.add(tokens[0]);
            sampleToTissue.put(tokens[0], tokens[1]);
        }
        fu.close();
        
        fileName = DIR_NAME + "GSE6919_GPL8300_121813.txt";
//        fileName = DIR_NAME + "GSE6919_GPL93_121813.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName);
        // Check the tissue source
        boolean isCheckingDone = false;
        // Do a filtering so that only primary cancers are kept
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            if (!isCheckingDone) {
                for (String sample : sampleToValue.keySet())
                    System.out.println(sample + "\t" + sampleToTissue.get(sample));
                isCheckingDone = true;
            }
            sampleToValue.keySet().retainAll(primary);
        }
        // Just a quick check
        Map<String, Double> lfng = geneToSampleToValue.get("LFNG");
        Map<String, Double> nkx31 = geneToSampleToValue.get("NKX3-1");
        Collections.sort(primary);
        System.out.println("\nCheck values");
        System.out.println("Sample\tLFNG\tNKX3_1");
        for (String sample : primary) {
            System.out.println(sample + "\t" + lfng.get(sample) + "\t" + nkx31.get(sample));
        }
        System.out.println("\nGSE6919 LFNG Gene Expression Correlation:");
        checkLFNGGeneCorrelations(geneToSampleToValue);
    }
    
    @Test
    public void checkMSKCCLFGNCorrelations() throws Exception {
        // Check with primary samples only
        List<String> primarySamples = new ArrayList<String>();
        String fileName = DIR_NAME + "MSKCC/MSKCC_PCa_Clinical_Annotation_LFNG.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[1].equals("PRIMARY"))
                primarySamples.add(tokens[0]);
        }
        fu.close();
        fileName = DIR_NAME + "MSKCC/MSKCC_PCa_mRNA_Zscores_Modified.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName,
                                                                                                            false);
        // Filter to primary samples only
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            sampleToValue.keySet().retainAll(primarySamples);
        }
        checkLFNGGeneCorrelations(geneToSampleToValue);
    }
    
    @Test
    public void outputTCGAGleasonScoreLFNGAndNKX31() throws Exception {
        String fileName = TCGA_DIR + "PRAD.uncv2.mRNAseq_RSEM_normalized_log2.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(fileName, false);
        Map<String, Double> lfngExp = new HashMap<String, Double>();
        Map<String, Double> nkx31Exp = new HashMap<String, Double>();
        for (String gene : geneToSampleToValue.keySet()) {
            if (gene.startsWith("LFNG")) {
                if (lfngExp.size() > 0)
                    throw new IllegalStateException("More than on LFNG gene!");
                lfngExp = filterTCGAExp(geneToSampleToValue.get(gene));
            }
            else if (gene.startsWith("NKX3-1")) {
                if (nkx31Exp.size() > 0)
                    throw new IllegalStateException("More than one NKX3-1 gene!");
                nkx31Exp = filterTCGAExp(geneToSampleToValue.get(gene));
            }
        }
        String clinFileName = "datasets/TCGA/PRAD/gdac.broadinstitute.org_PRAD.Clinical_Pick_Tier1.Level_4.2013111400.0.0/" + 
                              "PRAD.clin.merged.picked.transformed.txt";
        // Load sample to gleanson score
        Map<String, Double> sampleToGleanson = new HashMap<String, Double>();
        fu.setInput(clinFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Multiple gleason scores existing
            Double score = new Double(tokens[22]);
            sampleToGleanson.put(tokens[0].toUpperCase(), score);
        }
        fu.close();
//        String destFileName = DIR_NAME + "TCGAGleasonLFNGNKX3-1_122313.txt";
        String destFileName = DIR_NAME + "TCGAGleasonLFNGNKX3-1_Normal_122913.txt";
        outputLFNGNKX31Gleason(lfngExp, 
                               nkx31Exp, 
                               sampleToGleanson,
                               destFileName);
    }

    private void outputLFNGNKX31Gleason(Map<String, Double> lfngExp,
                                        Map<String, Double> nkx31Exp,
                                        Map<String, Double> sampleToGleanson,
                                        String destFileName) throws IOException {
        // We want to get the largest list even though there are some NA values
        List<String> sampleList = new ArrayList<String>(sampleToGleanson.keySet());
        if (sampleList.size() < lfngExp.size())
            sampleList = new ArrayList<String>(lfngExp.keySet());
        if (sampleList.size() < nkx31Exp.size())
            sampleList = new ArrayList<String>(nkx31Exp.keySet());
        Collections.sort(sampleList);
        fu.setOutput(destFileName);
        fu.printLine("Sample\ttype\tgleason_score\tLFNG\tNKX3-1");
        for (String sample : sampleList) {
            int index = sample.indexOf(".");
            String sampleId = sample.substring(0, index);
            String type = sample.substring(index + 1);
            Double gleason = sampleToGleanson.get(sampleId);
            if (type.equals("normal"))
                gleason = null; // Don't assign anything to gleason if it is a normal sample
            Double lfng = lfngExp.get(sample);
            Double nkx31 = nkx31Exp.get(sample);
            fu.printLine(sampleId + "\t" + 
                         type + "\t" +
                         (gleason == null ? "NA" : gleason) + "\t" + 
                         (lfng == null ? "NA" : lfng) + "\t" + 
                         (nkx31 == null ? "NA" : nkx31));
        }
        fu.close();
    }
    
    private Map<String, Double> filterTCGAExp(Map<String, Double> sampleToValue) {
        Map<String, Double> rtn = new HashMap<String, Double>();
        for (String sample : sampleToValue.keySet()) {
            String sampleId = sample.substring(0, 12);
            if (sample.endsWith("-11")) {
                rtn.put(sampleId + ".normal",
                        sampleToValue.get(sample));
            }
            else
                rtn.put(sampleId + ".primary",
                        sampleToValue.get(sample));
        }
        return rtn;
    }
    
    /**
     * Do a correlation analysis between LFNG and multiple others
     * @throws Exception
     */
    @Test
    public void checkTCGAGeneCorrelations() throws Exception {
        String fileName = TCGA_DIR + "PRAD.uncv2.mRNAseq_RSEM_normalized_log2.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(fileName, false);
        // Remove normal samples
        for (String gene : geneToSampleToValue.keySet()) {
            for (Iterator<String> it = geneToSampleToValue.get(gene).keySet().iterator(); it.hasNext();) {
                String sample = it.next();
                if (sample.endsWith("-11"))
                    it.remove();
            }
        }
        checkLFNGGeneCorrelations(geneToSampleToValue);
    }

    private void checkLFNGGeneCorrelations(Map<String, Map<String, Double>> geneToSampleToValue)
            throws MathException {
        // Get LFNG gene expression
        Map<String, Double> lfng = null;
        for (String gene : geneToSampleToValue.keySet()) {
            if (gene.startsWith("LFNG")) {
                lfng = geneToSampleToValue.get(gene);
                break;
            }
        }
        String[] checking = new String[] {"LFNG", "MFNG", "NOTCH", "DLL", "JAG", "HES", "HEY", "NKX"};
        boolean needToCheck = false;
        for (String gene : geneToSampleToValue.keySet()) {
           needToCheck = false;
           for (String tmp : checking) {
               if (gene.startsWith(tmp)) {
                   needToCheck = true;
                   break;
               }
           }
           if (!needToCheck)
               continue;
           Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
           Set<String> sharedSamples = InteractionUtilities.getShared(lfng.keySet(),
                                                                      sampleToValue.keySet());
           List<Double> values1 = new ArrayList<Double>();
           List<Double> values2 = new ArrayList<Double>();
           for (String sample : sharedSamples) {
               Double value1 = lfng.get(sample);
               Double value2 = sampleToValue.get(sample);
               values1.add(value1);
               values2.add(value2);
           }
           PearsonsCorrelation corTest = MathUtilities.constructPearsonCorrelation(values1, values2);
           System.out.println("Correlation between LFNG and " + gene + " with data points " + values1.size());
           System.out.println(corTest.getCorrelationMatrix().getEntry(0, 1) + 
                              " with pvalue " + corTest.getCorrelationPValues().getEntry(0, 1));
           System.out.println("----------------------------------------\n");
        }
    }
    
    @Test
    public void checkLFNGExpressionInTCGA() throws Exception {
        String fileName = DIR_NAME + "TCGA_LFNG_Exp.txt";
        fu.setInput(fileName);
        String line1 = fu.readLine();
        String[] headers = line1.split("\t");
        System.out.println("Total samples: " + (headers.length - 1));
        String line2 = fu.readLine();
        String[] values = line2.split("\t");
        System.out.println("Total values: " + (values.length - 1));
        List<Double> cancer = new ArrayList<Double>();
        List<Double> normal = new ArrayList<Double>();
        for (int i = 1; i < headers.length; i++) {
            if (headers[i].endsWith("-11"))
                normal.add(new Double(values[i]));
            else if (headers[i].endsWith("-01"))
                cancer.add(new Double(values[1]));
        }
        fu.close();
        System.out.println("Normal: " + normal.size());
        calculateMean(normal);
        System.out.println("Primary Cancer: " + cancer.size());
        calculateMean(cancer);
        // Do some t-test
        System.out.println("T-test between normal and primary:");
        double pvalue = MathUtilities.calculateTTest(normal, cancer);
        System.out.println("Pvalue: " + pvalue);
    }
    
    @Test
    public void checkLFNGExpressionInMSKCCDataSet() throws Exception {
        String fileName = DIR_NAME + "MSKCC_LFNG_Exp.txt";
        fu.setInput(fileName);
        String line1 = fu.readLine();
        String[] headers = line1.split("\t");
        String line2 = fu.readLine();
        String[] values = line2.split("\t");
        List<Double> cancer = new ArrayList<Double>();
        List<Double> mets = new ArrayList<Double>();
        List<Double> normal = new ArrayList<Double>();
        for (int i = 2; i < headers.length; i++) {
            if (headers[i].startsWith("PCA")) {
                String number = headers[i].substring(4);
                if (new Integer(number) < 182)
                    cancer.add(new Double(values[i]));
                else
                    mets.add(new Double(values[i]));
            }
            else
                normal.add(new Double(values[i]));
        }
        fu.close();
        System.out.println("Normal: " + normal.size());
        calculateMean(normal);
        System.out.println("Primary Cancer: " + cancer.size());
        calculateMean(cancer);
        System.out.println("Met: " + mets.size());
        calculateMean(mets);
        // Do some t-test
        System.out.println("T-test between normal and primary:");
        double pvalue = MathUtilities.calculateTTest(normal, cancer);
        System.out.println("Pvalue: " + pvalue);
        System.out.println("T-test between normal and mets: ");
        pvalue = MathUtilities.calculateTTest(normal, mets);
        System.out.println("Pvalue: " + pvalue);
    }
    
    private void calculateMean(List<Double> values) {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (Double value : values)
            stat.addValue(value);
        System.out.println("Mean: " + stat.getMean());
    }
    
    @Test
    public void calculatePairWiseCorrelations() throws Exception {
        String fileName = DIR_NAME + "GSE6919_GPL93_Checking_Genes_121813.txt";
//        String fileName = DIR_NAME + "GSE6919_GPL8300_Checking_Genes_121813.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<GeneData> data = new ArrayList<GeneData>();
        Set<String> genes = new HashSet<String>();
        Set<String> types = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            GeneData data1 = new GeneData();
            data1.gene = tokens[0];
            data1.sample = tokens[1];
            data1.type = tokens[2];
            data1.value = new Double(tokens[3]);
            data.add(data1);
            genes.add(data1.gene);
            types.add(data1.type);
        }
        fu.close();
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        System.out.println("Gene: " + geneList);
        for (String gene : geneList)
            System.out.print(gene + "\t");
        System.out.println();
        for (String gene : geneList)
            System.out.println(gene);
        for (String type : types) {
            if (!type.equals("primary"))
                continue;
            System.out.println(type);
            StringBuilder corrs = new StringBuilder();
            StringBuilder pvalues = new StringBuilder();
            for (int i = 0; i < geneList.size(); i++) {
                String gene1 = geneList.get(i);
                Map<String, Double> sampleToValue1 = getSampleToValue(data, gene1, type);
                for (int j = 0; j < geneList.size(); j++) {
                    if (i == j) {
                        corrs.append("\t");
                        pvalues.append("\t");
                        continue;
                    }
                    String gene2 = geneList.get(j);
                    Map<String, Double> sampleToValue2 = getSampleToValue(data, gene2, type);
                    PearsonsCorrelation corr = calculateCorrelation(sampleToValue1, sampleToValue2);
                    if (corr == null) {
                        corrs.append("\t");
                        pvalues.append("\t");
                        continue;
                    }
                    System.out.println(gene1 + " " + gene2);
                    System.out.println(corr.getCorrelationMatrix().getEntry(0, 1) + ", " + corr.getCorrelationPValues().getEntry(0, 1));
                    corrs.append(corr.getCorrelationMatrix().getEntry(0, 1)).append("\t");
                    pvalues.append(corr.getCorrelationPValues().getEntry(0, 1)).append("\t");
                    System.out.println();
                }
                corrs.append("\n");
                pvalues.append("\n");
            }
            System.out.println();
            System.out.println(corrs.toString());
            System.out.println(pvalues.toString());
        }
    }
    
    private PearsonsCorrelation calculateCorrelation(Map<String, Double> sampleToValue1,
                                                      Map<String, Double> sampleToValue2) {
        Set<String> samples = new HashSet<String>(sampleToValue1.keySet());
        samples.retainAll(sampleToValue2.keySet());
        List<Double> list1 = new ArrayList<Double>();
        List<Double> list2 = new ArrayList<Double>();
        for (String sample : samples) {
            Double value1 = sampleToValue1.get(sample);
            list1.add(value1);
            Double value2 = sampleToValue2.get(sample);
            list2.add(value2);
        }
        if (list1.size() == 0 || list2.size() == 0)
            return null;
        System.out.println("Size: " + list1.size());
        return MathUtilities.constructPearsonCorrelation(list1, list2);
    }
    
    private Map<String, Double> getSampleToValue(List<GeneData> data,
                                                 String gene,
                                                 String type) {
        Map<String, Double> sampleToValue = new HashMap<String, Double>();
        for (GeneData geneData : data) {
            if (geneData.gene.equals(gene) && geneData.type.equals(type)) {
                sampleToValue.put(geneData.sample, geneData.value);
            }
        }
        return sampleToValue;
    }
    
    @Test
    public void generateGeneToSampleToValue() throws IOException {
        List<String> checkingGenes = getCheckGenes();
        Map<String, String> sampleToType = loadSampleToTissueType();
        Map<String, Set<String>> sourceToSamples = loadSourceToSamples();
        
        //String fileName = DIR_NAME + "GSE6919_Global_Z_022213.txt";
//        String fileName = DIR_NAME + "GSE6919_GPL93_121813.txt";
        String fileName = DIR_NAME + "GSE6919_GPL8300_121813.txt";
        fu.setInput(fileName);
        String line = fu.readLine(); 
        String[] tokens = line.split("\t");
        List<String> sampleList = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++)
            sampleList.add(tokens[i]);
//        fu.setOutput(DIR_NAME + "GSE6919_Checking_Genes.txt");
//        fu.setOutput(DIR_NAME + "GSE6919_Checking_Genes_022713.txt");
//        fu.setOutput(DIR_NAME + "GSE6919_Checking_Genes_040813.txt");
//        fu.setOutput(DIR_NAME + "GSE6919_Checking_Genes_Sample_Avg_042913.txt");
//        fu.setOutput(DIR_NAME + "GSE6919_GPL93_Checking_Genes_121813.txt");
        fu.setOutput(DIR_NAME + "GSE6919_GPL8300_Checking_Genes_121813.txt");
        boolean needSampleAverage = false;
        fu.printLine("Gene\tSample\tType\tValue");
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (checkingGenes.contains(tokens[0])) {
                if (needSampleAverage) {
                    outputWithAverage(tokens, 
                                      sampleToType,
                                      sourceToSamples,
                                      sampleList);
                    continue;
                }
                for (int i = 1; i < tokens.length; i++) {
                    String token = tokens[i];
                    if (token.length() == 0)
                        continue;
                    String sample = sampleList.get(i - 1);
                    String type = sampleToType.get(sample);
                    fu.printLine(tokens[0] + "\t" + sample + "\t" + type + "\t" + tokens[i]);
                }
//                break;
            }
        }
        fu.close();
    }
    
    private void outputWithAverage(String[] tokens,
                                   Map<String, String> sampleToType,
                                   Map<String, Set<String>> sourceToSamples,
                                   List<String> sampleList) throws IOException {
        // Load all values
        Map<String, Double> sampleToValue = new HashMap<String, Double>();
        for (int i = 1; i < tokens.length; i++) {
            String token = tokens[i];
            if (token.length() == 0)
                continue;
            String sample = sampleList.get(i - 1);
            sampleToValue.put(sample, 
                              new Double(tokens[i]));
        }
        // Need to do average
        for (String source : sourceToSamples.keySet()) {
            Set<String> samples = sourceToSamples.get(source);
            double total = 0.0d;
            int count = 0;
            for (String sample : samples) {
                Double value = sampleToValue.get(sample);
                if (value == null)
                    continue;
                total += value;
                count ++;
            }
            if (count == 0)
                continue;
            fu.printLine(tokens[0] + "\t" + 
                         source + "\t" + 
                         getTypesForSamples(samples, sampleToType) + "\t" + 
                         total / count);
        }
    }
    
    /**
     * Load sample to tissue types from a pre-generated sample information file.
     * @return
     * @throws IOException
     */
    private Map<String, String> loadSampleToTissueType() throws IOException {
        return loadSampleToFetaure("Type");
    }
    
    private List<String> loadGSE6919PrimarySamples() throws IOException {
        Map<String, String> sampleToType = loadSampleToTissueType();
        List<String> rtn = new ArrayList<String>();
        for (String sample : sampleToType.keySet()) {
            String type = sampleToType.get(sample);
            if (type.equals("primary"))
                rtn.add(sample);
        }
        return rtn;
    }
    
    private Map<String, String> loadSampleToFetaure(String featureName) throws IOException {
        String fileName = DIR_NAME + "GSE6919_Samples_121913.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int index = -1;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals(featureName)) {
                index = i;
                break;
            }
        }
        if (index == -1)
            throw new IllegalArgumentException("Unknown feature name: " + featureName);
        Map<String, String> sampleToType = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            sampleToType.put(tokens[0], tokens[index]);
        }
        fu.close();
        return sampleToType;
    }
    
    /**
     * Generate sample information from a soft family file.
     * @throws IOException
     */
    @Test
    public void generateSampleInformation() throws IOException {
        String fileName = DIR_NAME + "GSE6919_family.soft";
        Map<String, String> sampleToTissue = dataHandler.loadSampleFeature(fileName,
                                                                           "!Sample_characteristics_ch1 = Tissue",
                                                                           ":");
        Map<String, String> sampleToSource = dataHandler.loadSampleFeature(fileName,
                                                                           "!Sample_source_name_ch1",
                                                                           "=");
        Map<String, String> sampleToGrade = dataHandler.loadSampleFeature(fileName,
                                                                          "!Sample_characteristics_ch1 = Gleason Grade", 
                                                                          ":");
        Map<String, String> sampleToStage = dataHandler.loadSampleFeature(fileName,
                                                                          "!Sample_characteristics_ch1 = Tumor stage",
                                                                          ":");
        System.out.println("Total samples: " + sampleToTissue.size());
        List<String> sampleList = new ArrayList<String>(sampleToTissue.keySet());
        Collections.sort(sampleList);
        System.out.println("Sample\tSource\tTissue\tType\tGleason Grade\tTumor stage");
        for (String sample : sampleList) {
            String source = sampleToSource.get(sample);
            String tissue = sampleToTissue.get(sample);
            String type = null;
            if (tissue.contains("metastases"))
                type = "metastases";
            else if (tissue.equals("primary prostate tumor"))
                type = "primary";
            else if (tissue.equals("normal prostate tissue free of any pathological alteration from brain-dead organ donor"))
                type = "normal_free";
            else if (tissue.equals("normal prostate tissue adjacent to tumor"))
                type = "normal_adjacent";
            String grade = sampleToGrade.get(sample);
            if (grade == null)
                grade = "NA";
            String stage = sampleToStage.get(sample);
            if (stage == null)
                stage = "NA";
            System.out.println(sample + "\t" + source + "\t" + tissue + "\t" + type + 
                               "\t" + grade + "\t" + stage);
        }
        Set<String> tissueTypes = new HashSet<String>(sampleToTissue.values());
        System.out.println("\nTissue types: " + tissueTypes.size());
        List<String> tissueTypeList = new ArrayList<String>(tissueTypes);
        Collections.sort(tissueTypeList);
        for (String type : tissueTypeList)
            System.out.println(type);
    }
    
    @Test
    //TODO: Merge same samples but duplicated (this may be done later on)
    public void processGSE6919SoftFamilyFile() throws IOException {
        dataHandler.setControlColName("SPOT_ID");
        dataHandler.setControlTypeName("--Control");
        
        
        String fileName = DIR_NAME + "GSE6919_family.soft";
//        String outFileName = DIR_NAME + "GSE6919_Global_Z_022213.txt";
        //String outFileName = DIR_NAME + "GSE6919_Global_Z_022713.txt";
//        String outFileName = DIR_NAME + "GSE6919_GPL93_121813.txt";
        String outFileName = DIR_NAME + "GSE6919_GPL8300_121813.txt";
        List<String> platforms = dataHandler.loadPlatforms(fileName);
        System.out.println("Total platforms: " + platforms.size() + " " + platforms);
        Map<String, String> probeIdToGene = dataHandler.loadProbeIdToGene(fileName);
        System.out.println("Total probe ids: " + probeIdToGene.size());
        platforms.clear();
//        platforms.add("GPL93");
        platforms.add("GPL8300");
        dataHandler.processSoftTextFile(fileName, 
                                        outFileName, 
                                        probeIdToGene, 
                                        platforms, 
                                        false);
    }
    
    @Test
    public void checkSampleInfo() throws IOException {
        String fileName = DIR_NAME + "GSE6919_Samples.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, String> sampleToSource = new HashMap<String, String>();
        Map<String, Set<String>> sourceToSamples = new HashMap<String, Set<String>>();
        Map<String, String> sampleToType = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int index = tokens[1].indexOf(",");
            String source = tokens[1].substring(index + 1).trim();
            // Need to get rid the last quotation mark
            source = source.substring(0, source.length() - 1);
            sampleToSource.put(tokens[0], source);
            InteractionUtilities.addElementToSet(sourceToSamples, 
                                                 source, 
                                                 tokens[0]);
            sampleToType.put(tokens[0],
                             tokens[3]);
        }
        fu.close();
        Set<String> samples = new HashSet<String>(sampleToSource.keySet());
        Set<String> sources = new HashSet<String>(sampleToSource.values());
        System.out.println("Total samples: " + samples.size());
        System.out.println("Total sources: " + sources.size());
        
        for (String source : sourceToSamples.keySet()) {
            samples = sourceToSamples.get(source);
            String type = getTypesForSamples(samples,
                                                   sampleToType);
            System.out.println(source + "\t" + 
                               samples.size() + "\t" +
                               samples + "\t" + 
                               type);
        }
    }
    
    private Map<String, Set<String>> loadSourceToSamples() throws IOException {
        String fileName = DIR_NAME + "GSE6919_Samples.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Set<String>> sourceToSamples = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            int index = tokens[1].indexOf(",");
            String source = tokens[1].substring(index + 1).trim();
            // Need to get rid the last quotation mark
            source = source.substring(0, source.length() - 1);
            InteractionUtilities.addElementToSet(sourceToSamples, 
                                                 source, 
                                                 tokens[0]);
        }
        fu.close();
        return sourceToSamples;
    }
    
    private String getTypesForSamples(Set<String> samples,
                                           Map<String, String> sampleToType) {
        Set<String> types = new HashSet<String>();
        for (String sample : samples)
            types.add(sampleToType.get(sample));
        if (types.size() > 1)
            throw new IllegalStateException(samples + " have more than two types: " + types);
        return types.iterator().next();
    }
    
    private List<String> getCheckGenes() {
        // DLL1 and NOTCH4 are not in GSE6919 data set.
        String[] genes = new String[] {
                "LFNG", "MFNG", "NOTCH1", "NOTCH2", "NOTCH3", "NOCTH4", "DLL1", "JAG1", "HES1", "HEY1", "HEY2", "HEYL", "NKX3-1"
        };
        return Arrays.asList(genes);
    }
    
    private class GeneData {
        String gene;
        String sample;
        String type;
        double value;
    }
}
