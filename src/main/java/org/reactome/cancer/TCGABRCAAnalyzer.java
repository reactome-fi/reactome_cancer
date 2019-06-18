/*
 * Created on Sep 2, 2014
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.junit.Test;
import org.reactome.pathway.factorgraph.GibbsSampler;
import org.reactome.pathway.factorgraph.ReactomePathwayFGRunner;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.ReactomeDBBridge;

/**
 * Perform some analyses related to TCGA BRCA data set.
 * @author gwu
 *
 */
public class TCGABRCAAnalyzer {
//    private final String DIR_NAME = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/paradigm/twoCases/brca/Node11/";
//    private final String DIR_NAME = "/Users/gwu/Documents/EclipseWorkspace/caBigR3/results/paradigm/twoCases/brca/102414/";
    // HNSC results
    // Only gene expression data is used in this run
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/121914_1/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/121014/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/121814/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/112914/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/112614/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/111714/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/111414/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/110614/";
//    private final String DIR_NAME = "results/paradigm/twoCases/hnsc/020215/";
    // Some test
    private final String DIR_NAME = "tmp/hnsc_gaussian/";
    
    public TCGABRCAAnalyzer() {
    }
    
    /**
     * Just a quick way to check correlation between two genes expression.
     * @throws Exception
     */
    @Test
    public void checkGeneExpressionCorrelations() throws Exception {
//        String fileName = "test_data/tcga_brca/brca.mRNA.txt";
        // Repeat a new file
        String fileName = "datasets/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014120600.0.0/" +
                "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.transformed.txt";

        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(fileName, false);
        String[] checkGenes = new String[] {
                "MFNG",
                "PIK3CG",
                "NOTCH4"
        };
        for (int j = 0; j < checkGenes.length; j++) {
            String gene1 = checkGenes[j];
            Map<String, Double> sampleToValue1 = geneToSampleToValue.get(gene1);
            for (int i = j + 1; i < checkGenes.length; i++) {
                String gene2 = checkGenes[i];
                Map<String, Double> sampleToValue2 = geneToSampleToValue.get(gene2);
                if (sampleToValue2 == null) {
                    System.out.println("Cannot find expression for " + gene2);
                    continue;
                }
                List<Double> valueList1 = new ArrayList<Double>();
                List<Double> valueList2 = new ArrayList<Double>();
                Set<String> sharedSamples = new HashSet<String>(sampleToValue1.keySet());
                sharedSamples.retainAll(sampleToValue2.keySet());
                for (String sample : sharedSamples) {
                    Double value1 = sampleToValue1.get(sample);
                    valueList1.add(value1);
                    Double value2 = sampleToValue2.get(sample);
                    valueList2.add(value2);
                }
                PearsonsCorrelation cor = MathUtilities.constructPearsonCorrelation(valueList1, valueList2);
                System.out.println(gene1 + "\t" + gene2);
                System.out.println("Cor: " + cor.getCorrelationMatrix().getEntry(0, 1));
                System.out.println("P-value: " + cor.getCorrelationPValues().getEntry(0, 1));
                System.out.println("Total shared samples: " + sharedSamples.size());
                System.out.println();
            }
        }
    }
    
    @Test
    public void checkSharedGenesInTopPathways() throws Exception {
//        List<GKInstance> allPathways = new PARADIGMRunner().getPathwayList();
//        PersistenceAdaptor dba = allPathways.get(0).getDbAdaptor();
//        String fileName = DIR_NAME + "PathwaySummary_OSTypes_NoMeta_All_Combined_FDR.txt";
//        // Get the top pathways with FDRs < 0.05
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = fu.readLine();
//        List<GKInstance> pathways = new ArrayList<GKInstance>();
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            Double fdr = new Double(tokens[tokens.length - 1]);
//            if (fdr > 0.05d)
//                break;
//            Long dbId = new Long(tokens[0]);
//            GKInstance pathway = dba.fetchInstance(dbId);
//            pathways.add(pathway);
//        }
//        fu.close();
//        Set<String> sharedGenes = new HashSet<String>();
//        System.out.println("\nShared genes between pathways having FDR < 0.05:");
//        for (int i = 0; i < pathways.size() - 1; i++) {
//            GKInstance pathway1 = pathways.get(i);
//            Set<String> genes1 = getGenesInPathway(pathway1);
//            for (int j = i + 1; j < pathways.size(); j++) {
//                GKInstance pathway2 = pathways.get(j);
//                Set<String> genes2 = getGenesInPathway(pathway2);
//                Set<String> shared = InteractionUtilities.getShared(genes1, genes2);
//                System.out.println(pathway1.getDisplayName() + "\t" + 
//                                   pathway2.getDisplayName() + "\t" + 
//                                   shared.size() + "\t" + 
//                                   shared);
//                sharedGenes.addAll(shared);
//            }
//        }
//        System.out.println("\nShared genes in other pathways:");
//        for (GKInstance pathway : allPathways) {
//            Set<String> genes = getGenesInPathway(pathway);
//            Set<String> shared = InteractionUtilities.getShared(genes, sharedGenes);
//            if (shared.size() == 0)
//                continue;
//            System.out.println(pathway + ": " + shared);
//        }
    }
    
    private Set<String> getGenesInPathway(GKInstance pathway) throws Exception {
        Set<GKInstance> refSeqs = InstanceUtilities.grepRefPepSeqsFromPathway(pathway);
        Set<String> genes = new HashSet<String>();
        for (GKInstance refSeq : refSeqs) {
            if (!refSeq.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
                continue;
            }
            String gene = (String) refSeq.getAttributeValue(ReactomeJavaConstants.geneName);
            if (gene != null)
                genes.add(gene);
        }
        return genes;
    }
    
    private List<File> getIPAFiles() throws IOException {
        String dirName = DIR_NAME;
        File dir = new File(dirName);
        List<File> ipaFiles = new ArrayList<File>();
        
        for (File file : dir.listFiles()) {
            String fileName = file.getName();
            if (file.isDirectory() && fileName.startsWith("Node")) {
                for (File file1 : file.listFiles()) {
                    fileName = file1.getName();
                    if (fileName.endsWith(".txt") && 
                        !fileName.equals("logging.txt") && 
                        !fileName.equals("Pathways.txt"))
                        ipaFiles.add(file1);
                }
            }
        }
        
//        // Some test files
//        String[] fileNames = new String[] {
////                "IPA Transcriptional_Regulation_of_Pluripotent_Stem_Cells.txt",
////                "new/IPA Transcriptional_Regulation_of_Pluripotent_Stem_Cells.txt"
////                "Transcriptional_Regulation_of_Pluripotent_Stem_Cells.txt"
////                "Node8/Passive_transport_by_Aquaporins.txt"
////                "YAP1-_and_WWTR1__TAZ_-stimulated_gene_expression.txt"
////                "Transcriptional_regulation_of_pluripotent_stem_cells.txt"
////                "Metabolism_of_nucleotides.txt"
//                "Transcriptional_regulation_of_pluripotent_stem_cells.txt"
//        };
//        for (String fileName : fileNames) {
//            File file = new File(dirName, fileName);
//            ipaFiles.add(file);
//        }
        
        return ipaFiles;
    }
    
    private List<String> getFinishedPathways() throws IOException {
        File file = new File(DIR_NAME, "RunningTimes.txt");
        if (!file.exists())
            return null;
        FileUtility fu = new FileUtility();
        fu.setInput(file.getAbsolutePath());
        String line = null;
        List<String> rtn = new ArrayList<String>();
        boolean inTheList = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("Pathways running time (minutes):")) {
                inTheList = true;
            }
            else if (inTheList) {
                if (line.startsWith("In total"))
                    break;
                int index = line.indexOf("\t");
                String pathway = line.substring(0, index);
                rtn.add(pathway);
            }
        }
        fu.close();
        return rtn;
    }
    
    @Test
    public void runGibbsLearningAndAnalyzerResults() throws Exception {
        GibbsSampler sampler = new GibbsSampler();
        sampler.learnParametersFromGeneratedSamples();
        analyzeFactorGraphResults();
    }
    
    /**
     * This method is used to check the analysis results from the PARADIGM method without running the background.
     * @throws IOException
     */
    @Test
    public void analyzeFactorGraphResults() throws Exception {
        FileUtility.initializeLogging();
        String dirName = DIR_NAME;
        List<File> ipaFiles = getIPAFiles();
        System.out.println("Total IPA files: " + ipaFiles.size());
        // Get the type information
//        Map<String, Set<String>> typeToSamples = loadTypeToSamples();
//        String summaryFileName = dirName + "PathwaySummary.txt";
//        String detailedInfoFileName = dirName + "PathwayDetailedInfo.txt";
        
        // Use type information based on overall survival for dead patients only
//        Map<String, Set<String>> typeToSamples = loadOSTypeToSamples(false);
////        String summaryFileName = dirName + "PathwaySummary_OSTypes.txt";
////        String detailedInfoFileName = dirName + "PathwayDetailedInfo_OSTypes.txt";
//        String summaryFileName = dirName + "PathwaySummary_OSTypes_All_Combined.txt";
//        String detailedInfoFileName = dirName + "PathwayDetailedInfo_OSTypes_All_Combined.txt";
        // As above, but excluding metastatic samples
//        Map<String, Set<String>> typeToSamples = loadOSTypeToSamples(true);
//        String summaryFileName = dirName + "PathwaySummary_OSTypes_NoMeta.txt";
//        String detailedInfoFileName = dirName + "PathwayDetailedInfo_OSTypes_NoMeta.txt";
        // Test with combined p-values from all outputs
//        String summaryFileName = dirName + "PathwaySummary_OSTypes_NoMeta_All_Combined.txt";
//        String detailedInfoFileName = dirName + "PathwayDetailedInfo_OSTypes_NoMeta_All_Combined.txt";
        
        // Results for HNSC
        String summaryFileName = dirName + "PathwaySummary_Progressor_All_Combined.txt";
        String detailedInfoFileName = dirName + "PathwayDetailedInfo_Progressor_All_Combined.txt";
        Map<String, Set<String>> typeToSamples = loadHNSCTypeToSamples();
        
//        // For synthetic data
//        Map<String, Set<String>> typeToSamples = loadGeneratedSamples(false);
        
        // Just a quick check
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            System.out.println(type + ": " + samples.size());
        }
        if (true)
            return;
        // Used to map from file name to pathway
        Map<String, GKInstance> fileNameToPathway = getMapFromFileNameToPathway();
        Map<String, List<Double>> idToValues0 = new HashMap<String, List<Double>>();
        Map<String, List<Double>> idToValues1 = new HashMap<String, List<Double>>();
        FileUtility summaryFu = new FileUtility();
        summaryFu.setOutput(summaryFileName);
        summaryFu.printLine("DBID\tPathwayName\tUp\tDown\tCombinedPValue\tMinP");
        // Hold the detailed information for a pathway
        FileUtility detailedFu = new FileUtility();
        detailedFu.setOutput(detailedInfoFileName);
        detailedFu.printLine("DBID\tPathwayName\tOutputId\tOutputName\tLowerMean\tUpperMean\tDiff\tP-value");
        List<String> finishedPathways = getFinishedPathways();
        System.out.println();
        for (String fileName : fileNameToPathway.keySet()) {
            System.out.println(fileName + "\t" + fileNameToPathway.get(fileName));
        }
        System.out.println();
        for (File file : ipaFiles) {
            String fileName = file.getName();
            // For some reason, Mac OS change the file names
            GKInstance pathway = fileNameToPathway.get(fileName);
            if (pathway == null)
                pathway = fileNameToPathway.get(fileName.toLowerCase());
            if (pathway == null)
                throw new IllegalStateException(file.getAbsolutePath() + " cannot be mapped to a pathway!");
            if (finishedPathways != null && !finishedPathways.contains(pathway.toString())) {
                System.out.println(pathway + " has not finished yet!");
                continue;
            }
            idToValues0.clear();
            idToValues1.clear();
            if(!loadIPAValues(file, 
                             pathway,
                             typeToSamples, 
                             idToValues0,
                             idToValues1))
                continue;
            outputIPAComparison(pathway,
                                idToValues0,
                                idToValues1,
                                detailedFu,
                                summaryFu);
        }
        summaryFu.close();
        detailedFu.close();
        // Calculate FDR values
        int index = summaryFileName.lastIndexOf(".");
        String fdrFileName = summaryFileName.substring(0, index) + "_FDR" + summaryFileName.substring(index);
        calculateFDRs(summaryFileName, fdrFileName);
    }
    
    private void calculateFDRs(String srcFileName,
                               String targetFileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        String header = fu.readLine();
        String[] tokens = header.split("\t");
        List<String> lines = new ArrayList<String>();
        String line = null;
        while ((line = fu.readLine()) != null)
            lines.add(line);
        fu.close();
        Map<String, Double> combinedPFDRs = calculateFDRs(lines, 
                                                          tokens.length  - 2);
        Map<String, Double> minPFDRs = calculateFDRs(lines,
                                                     tokens.length - 1);
        // Output
        fu.setOutput(targetFileName);
        fu.printLine(header + "\tFDR_CombinedP\tFDR_minP");
        for (int i = 0; i < lines.size(); i++) {
            fu.printLine(lines.get(i) + "\t" + 
                         combinedPFDRs.get(lines.get(i)) + "\t" +
                         minPFDRs.get(lines.get(i)));
        }
        fu.close();
        // Delete the original file.
        new File(srcFileName).delete();
    }

    private Map<String, Double> calculateFDRs(List<String> lines,
                                              final int pvalueCol) {
        // Do a sorting based on p-values
        Collections.sort(lines, new Comparator<String>() {
            public int compare(String line1, String line2) {
                String[] tokens = line1.split("\t");
                Double pvalue1 = new Double(tokens[pvalueCol]);
                tokens = line2.split("\t");
                Double pvalue2 = new Double(tokens[pvalueCol]);
                return pvalue1.compareTo(pvalue2);
            }
        });
        List<Double> pvalues = new ArrayList<Double>();
        for (String line1 : lines) {
            String[] tokens = line1.split("\t");
            Double pvalue = new Double(tokens[pvalueCol]);
            pvalues.add(pvalue);
        }
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        Map<String, Double> lineToPvalue = new HashMap<String, Double>();
        for (int i = 0; i < lines.size(); i++) {
            lineToPvalue.put(lines.get(i), fdrs.get(i));
        }
        return lineToPvalue;
    }
    
    private void outputIPAComparison(GKInstance pathway,
                                     Map<String, List<Double>> idToValues0,
                                     Map<String, List<Double>> idToValues1,
                                     FileUtility detailedFu,
                                     FileUtility summaryFu) throws Exception {
        MannWhitneyUTest uTest = new MannWhitneyUTest();
        List<Double> pvalues = new ArrayList<Double>();
        int upCount = 0;
        int downCount = 0;
        double upSum = 0.0d;
        double downSum = 0.0d;
        for (String key : idToValues0.keySet()) {
            List<Double> values0 = idToValues0.get(key);
            List<Double> values1 = idToValues1.get(key);
//            if (values0 == null || values1 == null || values0.size() == 0 || values1.size() == 0)
//                continue; // For some reason!
            double mean0 = MathUtilities.calculateMean(values0);
            double mean1 = MathUtilities.calculateMean(values1);
            // If there is an infinity value, just escape it
            if (Double.isInfinite(mean0) || Double.isInfinite(mean1))
                continue;
            double diff = mean1 - mean0;
            double pvalue = uTest.mannWhitneyUTest(convertDoubleListToArray(values1),
                                                   convertDoubleListToArray(values0));
            GKInstance output = pathway.getDbAdaptor().fetchInstance(new Long(key));
            if (output == null)
                throw new IllegalStateException("Cannot find instance for " + key);
            detailedFu.printLine(pathway.getDBID() + "\t" + 
                    pathway.getDisplayName() + "\t" + 
                    key + "\t" + 
                    output.getDisplayName() + "\t" + 
                    mean0 + "\t" + 
                    mean1 + "\t" + 
                    diff + "\t" + 
                    pvalue);
            // Combined p-values without filtering
            if (pvalue < 1.0e-14d)
                pvalue = 1.0e-14d;
            pvalues.add(pvalue);
            if (diff > 0.0d) {
                upSum += diff;
                upCount ++;
            }
            else if (diff < 0.0d) {
                downSum += diff;
                downCount ++;
            }
//            if (pvalue > 0.01d || Math.abs(diff) < 0.1d)
//                continue;
//            if (diff > 0.0d)
//                upCount ++;
//            else
//                downCount ++;
//            if (pvalue < 1.0e-14d)
//                pvalue = 1.0e-14d;
//            pvalues.add(pvalue);
        }
        //        double combinedPValue = new PValueCombiner().combinePValue(new ArrayList<List<Double>>(idToValues1.values()),
        //                                                                   pvalues);
        Double combinedPValue = null;
        if (pvalues.size() == 0)
            combinedPValue = 1.0d; // Use the maximum values for FDR calculation
        else
            combinedPValue = MathUtilities.combinePValuesWithFisherMethod(pvalues);
        Double minP = null;
        if (pvalues.size() == 0)
            minP = 1.0d;
        else
            minP = getMinimalPValue(pvalues);
        String text = (pathway.getDBID() + "\t" + 
                pathway.getDisplayName() + "\t" + 
                (upCount == 0 ? "0.0" : upSum / upCount) + "\t" + 
                (downCount == 0 ? "0.0" : downSum / downCount) + "\t" + 
                combinedPValue + "\t" + 
                minP);
        summaryFu.printLine(text);
    }
    
    private double getMinimalPValue(List<Double> pvalues) {
        double minP = Double.MAX_VALUE;
        for (Double pvalue : pvalues) {
            minP = Math.min(pvalue, minP);
        }
        return minP;
    }
    
    /**
     * Results here are compared with results generated by a R package, MADAM, for calculating
     * combined p-values from Fisher's method.
     * @throws MathException
     */
    @Test
    public void testCombinedPValue() throws MathException {
        String text = "0.4089769 0.551435 0.1029247 0.001";
        text = "0.0455565 0.4533342 0.04205953 0.0010000";
        List<Double> pvalues = new ArrayList<Double>();
        String[] tokens = text.split(" ");
        for (String token : tokens)
            pvalues.add(new Double(token));
        Double pvalue = MathUtilities.combinePValuesWithFisherMethod(pvalues);
        System.out.println("Combined: " + pvalue);
    }
    
    private double[] convertDoubleListToArray(List<Double> list) {
        double[] rtn = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            Double value = list.get(i);
            if (value == null)
                throw new IllegalArgumentException("Double List contains a null value!");
            rtn[i] = value;
        }
        return rtn;
    }

    private boolean loadIPAValues(File file,
                                  GKInstance pathway,
                                  Map<String, Set<String>> typeToSamples,
                                  Map<String, List<Double>> idToValues0,
                                  Map<String, List<Double>> idToValues1) throws Exception {
        // Check if there is anything in the file
        if (file.length() == 0) {
            System.err.println(file + " is empty!");
            return false; 
        }
        ReactomeDBBridge reactomeHelper = new ReactomeDBBridge();
        CancerGeneExpressionCommon dataHelper = new CancerGeneExpressionCommon();
        Set<Long> outputIds = reactomeHelper.getAllOutputIds(pathway);
//        System.out.println("File: " + file.getAbsolutePath());
        Map<String, Map<String, Double>> sampleToIdToValue = dataHelper.loadGeneExp(file.getAbsolutePath(), false);
        // Split values into two groups
        int type0 = 0, type1 = 0;
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            Map<String, List<Double>> idToValues = null;
            if (type.equals("0"))
                idToValues = idToValues0;
            else
                idToValues = idToValues1;
            for (String sample : samples) {
                Map<String, Double> idToValue = sampleToIdToValue.get(sample);
                if (idToValue == null) {
//                    System.err.println(sample + " has no value!");
                    continue; // It may not have value
                }
                if (type.equals("0"))
                    type0 ++;
                else
                    type1 ++;
                for (Long id : outputIds) {
                    Double value = idToValue.get(id.toString());
                    if (value == null)
                        continue;
                    List<Double> values = idToValues.get(id.toString());
                    if (values == null) {
                        values = new ArrayList<Double>();
                        idToValues.put(id.toString(), values);
                    }
                    values.add(value);
                }
            }
        }
        System.out.println(file + ": type 0 samples, " + type0 + "; type 1 samples " + type1 + "; ratio " + (double) type0 / (type0 + type1));
        return true;
    }
    
    private Map<String, GKInstance> getMapFromFileNameToPathway() throws Exception {
        ReactomePathwayFGRunner runner = new ReactomePathwayFGRunner();
//        MySQLAdaptor dba = new MySQLAdaptor("localhost", "reactome_51_plus_i", "root", "macmysql01");
//        runner.setAdaptor(dba);
        List<GKInstance> pathways = runner.getPathwayList();
        Map<String, GKInstance> rtn = new HashMap<String, GKInstance>();
        for (GKInstance pathway : pathways) {
            System.out.println(pathway);
            String pathwayFileName = InteractionUtilities.getFileNameFromInstance(pathway) + ".txt";
            System.out.println(pathwayFileName);
            String fileName = "IPA " + pathwayFileName;
            // Two cases
            rtn.put(pathwayFileName, pathway);
            rtn.put(fileName, pathway);
            rtn.put(pathwayFileName.toLowerCase(), pathway); // For some weird behavior under mac: file name cases are changed.
        }
        return rtn;
    }
    
    /**
     * Breast cancer gene signature created by using the MCL method based on microarray data sets.
     * @return
     */
    private List<String> getSignatureGenes() {
        // All 31 genes
        String text = "BIRC5,AURKB,KIF20A,NCAPD2,CENPN,BUB1,SPC25,RCC2,CDCA8,MAD2L1,NSUN2,SEC13,CCDC99,"
                + "NDC80,BUB1B,ANKZF1,KNTC1,ZWILCH,"
                + "CENPE,ZW10,CENPQ,MAD1L1,ERCC6L,AURKC,"
                + "CENPH,NUDC,INCENP,NCAPG,CENPI,TAOK1,ITGB3BP";
//        // The smaller signature contains 7 genes only
//        text = "ARUKB,BIRC5,BUB1,CENPN,KIF20A,NCAPD2,SPC25";
        List<String> list = new ArrayList<String>();
        String[] tokens = text.split(",");
        for (String token : tokens)
            list.add(token);
        return list;
    }
    
    /**
     * This method is used to generate files needed by GSEA.
     * @throws IOException
     */
    @Test
    public void generateFilesForGSEA() throws Exception {
        Map<String, Set<String>> typeToSamples = loadOSTypeToSamples(true);
        // Create a map from sample to type
        Map<String, String> sampleToType = new HashMap<String, String>();
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            for (String sample : samples)
                sampleToType.put(sample, type);
        }
        // Gene expression data
        String fileName = "test_data/tcga_brca/brca.mRNA.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // Gene expression file
        String outfileName = DIR_NAME + "brca.mRNA.GSEA.txt";
        fu.setOutput(outfileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        StringBuilder builder = new StringBuilder();
        builder.append("NAME\tDESCRIPTION");
        for (int i = 1; i < tokens.length; i++) {
            if (sampleToType.keySet().contains(tokens[i]))
                builder.append("\t").append(tokens[i]);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        String[] samples = tokens;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            builder.append(tokens[0]);
            builder.append("\tna"); // Nothing in description
            for (int i = 1; i < tokens.length; i++) {
                if (sampleToType.keySet().contains(samples[i])) {
                    String value = tokens[i];
                    if (value.equals("NA"))
                        value = ""; // required by GSEA
                    builder.append("\t").append(value);
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
        // Generate the class type
        outfileName = DIR_NAME + "brca.GSEA.cls";
        fu.setOutput(outfileName);
        int count = 0;
        for (int i = 1; i < samples.length; i++) {
            String type = sampleToType.get(samples[i]);
            if (type == null)
                continue;
            count ++;
            builder.append(type).append(" ");
        }
        fu.printLine(count + " 2 1");
        fu.printLine("# LOW UP");
        fu.printLine(builder.toString());
        builder.setLength(0);
        fu.close();
        // Generate GeneSet file
//        List<GKInstance> pathways = new PARADIGMRunner().getPathwayList();
//        outfileName = DIR_NAME + "reactome.GSEA.gmt";
//        fu.setOutput(outfileName);
//        for (GKInstance pathway : pathways) {
//            Set<String> genes = getGenesInPathway(pathway);
//            builder.append(pathway.getDisplayName() + "\tna");
//            for (String gene : genes)
//                builder.append("\t").append(gene);
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
    }
    
    @Test
    public void splitSamples() throws IOException {
        List<String> geneList = getSignatureGenes();
        System.out.println("Total genes in the signature: " + geneList.size());
        // Load the Firehose breast cancer RNA seq data
        String fileName = "test_data/tcga_brca/brca.mRNA.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = new CancerGeneExpressionCommon().loadGeneExp(fileName, false);
        geneToSampleToValue.keySet().retainAll(geneList); // Just focus on these genes in the signature
        System.out.println("Filter the dataset to the signature: " + geneToSampleToValue.size());
        // Get all samples first
        Set<String> samples = new HashSet<String>();
        for (String gene : geneToSampleToValue.keySet()) {
            samples.addAll(geneToSampleToValue.get(gene).keySet());
        }
        System.out.println("Total samples related to the signature: " + samples.size());
        SummaryStatistics stats = new SummaryStatistics();
        Map<String, Double> sampleToScore = new HashMap<String, Double>();
        for (String sample : samples) {
            Double score = calculateSignatureScore(sample, geneToSampleToValue, geneList);
            System.out.println(sample + "\t" + score);
            sampleToScore.put(sample, score);
            stats.addValue(score);
        }
        System.out.println("Sample\tMCLSignatureScore\tType");
        double mean = stats.getMean();
        for (String sample : sampleToScore.keySet()) {
            Double score = sampleToScore.get(sample);
            System.out.println(sample + "\t" + 
                               score + "\t" + 
                               (score >= mean ? 1 : 0));
        }
    }
    
    @Test
    public void testLoadOSSurvialTypeToSamples() throws IOException {
        Map<String, Set<String>> typeToSamples = loadOSTypeToSamples(true);
        System.out.println("Number of types: " + typeToSamples.size());
        for (String type : typeToSamples.keySet())
            System.out.println(type + ": " + typeToSamples.get(type).size());
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            String label = null;
            if (type.equals("0"))
                label = "Shorter";
            else
                label = "Longer";
            for (String sample : samples)
                System.out.println(sample + "\t" + label);
        }
    }
    
    /**
     * Use this method to split BRCA samples into two groups based on the medium of
     * overall survival for patients with vital-status 1, which are dead. 
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadOSTypeToSamples(boolean excludeMetaSamples) throws IOException {
        String clinFileName = "datasets/TCGA/BRCA/clinical/gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2014071500.0.0/" + 
                              "BRCA.clin.merged.picked.transformed.stages.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(clinFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        List<String> list = new ArrayList<String>();
        for (String token : tokens)
            list.add(token);
        int osEventIndex = list.indexOf("OSEVENT");
        int osTimeIndex = list.indexOf("OSDURATION");
        int metaIndex = list.indexOf("pathology.M.stage");
        Map<String, Double> sampleToTime = new HashMap<String, Double>();
        DescriptiveStatistics stat = new DescriptiveStatistics();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (excludeMetaSamples && !tokens[metaIndex].equals("m0"))
                continue; // Exclude samples that are not m0, not meta detected
            if (!tokens[osEventIndex].equals("1"))
                continue; // Not dead patient or null valued patient
            Double value = new Double(tokens[osTimeIndex]);
            sampleToTime.put(tokens[0], value);
            stat.addValue(value);
        }
        fu.close();
        // Split samples based on the medium of OS times
        Map<String, Set<String>> osTypeToSamples = new HashMap<String, Set<String>>();
        double medium = stat.getPercentile(50.0d);
        for (String sample : sampleToTime.keySet()) {
            Double time = sampleToTime.get(sample);
            String type = (time < medium ? "0" : "1");
            InteractionUtilities.addElementToSet(osTypeToSamples, type, sample);
        }
        return osTypeToSamples;
    }
    
    /**
     * Load a pre-generated sample's type information file.
     */
    private Map<String, Set<String>> loadTypeToSamples() throws IOException {
        Map<String, Set<String>> typeToSamples = new HashMap<String, Set<String>>();
        FileUtility fu = new FileUtility();
        String fileName = "results/paradigm/twoCases/brca/SampleToMCLSignatureScore.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            InteractionUtilities.addElementToSet(typeToSamples, tokens[2], tokens[0]);
        }
        fu.close();
        return typeToSamples;
    }
    
    private Map<String, Set<String>> loadHNSCTypeToSamples() throws IOException {
        String fileName = "/Users/gwu/Dropbox/OHSU/HNSC_Project/hnsc_followup_annotation.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Set<String>> typeToSamples = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            InteractionUtilities.addElementToSet(typeToSamples, tokens[1], tokens[0]);
        }
        // We want to use 0 and 1 for types
        Map<String, Set<String>> rtn = new HashMap<String, Set<String>>();
        Set<String> samples0 = typeToSamples.get("NonProgressor");
        rtn.put("0", samples0);
        System.out.println("0 for NonProgressor");
        Set<String> samples1 = typeToSamples.get("Progressor");
        rtn.put("1", samples1);
        System.out.println("1 for Progressor");
        return rtn;
    }
    
    @Test
    public void processObservationDataForHNSC() throws Exception {
        Map<String, Set<String>> typeToSamples = loadHNSCTypeToSamples();
        for (String type : typeToSamples.keySet()) {
            System.out.println(type + ": " + typeToSamples.get(type));
        }
//        String srcFileName = "tmp/Metabolism of nucleotides_GeneExp.txt";
//        String targetFileName = "tmp/Metabolism of nucleotides_GeneExp_Annotated.txt";
//        String srcFileName = "tmp/Metabolism of nucleotides_CNV.txt";
//        String targetFileName = "tmp/Metabolism of nucleotides_CNV_Annotated.txt";
        String srcFileName = "tmp/Metabolism_of_nucleotides.txt";
        String targetFileName = "tmp/Metabolism_of_nucleotides_Annotated.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        String header = line;
        Map<String, String> sampleToLine = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            int index = line.indexOf("\t");
            String sample = line.substring(0, index);
            sampleToLine.put(sample, line);
        }
        fu.close();
        fu.setOutput(targetFileName);
        fu.printLine(header);
        Map<String, List<Double>> sampleToValues = new HashMap<String, List<Double>>();
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            for (String sample : samples) {
                List<Double> values = new ArrayList<Double>();
                line = sampleToLine.get(sample);
                fu.printLine(line);
                String[] tokens = line.split("\t");
                for (int i = 1; i < tokens.length; i++) {
                    if (tokens[i].equals("na"))
                        values.add(null);
                    else
                        values.add(new Double(tokens[i]));
                }
                sampleToValues.put(sample, values);
            }
        }
        StringBuilder builder = new StringBuilder();
        builder.append("Diff");
        List<Double> list1 = new ArrayList<Double>();
        List<Double> list0 = new ArrayList<Double>();
        StringBuilder pvalueText = new StringBuilder();
        pvalueText.append("PValue");
        for (int i = 1; i < header.split("\t").length; i++) {
            list1.clear();
            list0.clear();
            for (String type : typeToSamples.keySet()) {
                Set<String> samples = typeToSamples.get(type);
//                System.out.println(type + ": " + samples.size());
                for (String sample : samples) {
                    List<Double> values = sampleToValues.get(sample);
                    Double value = values.get(i - 1);
                    if (value == null)
                        continue;
                    if (type.equals("1"))
                        list1.add(value);
                    else
                        list0.add(value);
                }
            }
            double mean1 = MathUtilities.calculateMean(list1);
            double mean0 = MathUtilities.calculateMean(list0);
            double diff = mean1 - mean0;
            double pvalue = MathUtilities.calculateTTest(list1, list0);
            builder.append("\t").append(diff);
            pvalueText.append("\t").append(pvalue);
        }
        fu.printLine(builder.toString());
        fu.printLine(pvalueText.toString());
        fu.close();
    }
    
    private Map<String, Set<String>> loadGeneratedSamples(boolean needRandomSamples) throws IOException {
        return new GibbsSampler().getTypeToSamples(needRandomSamples);
    }
    
    @Test
    public void testLoadHNSCTypeToSamples() throws IOException {
        Map<String, Set<String>> typeToSamples = loadHNSCTypeToSamples();
        Map<String, String> typeToName = new HashMap<String, String>();
        typeToName.put("0", "NonProgressor");
        typeToName.put("1", "Progressor");
        for (String type : typeToSamples.keySet()) {
            Set<String> samples = typeToSamples.get(type);
            System.out.println(type + ": " + samples.size());
            for (String sample : samples)
                System.out.println(sample + "\t" + typeToName.get(type));
        }
    }
    
    private double calculateSignatureScore(String sample,
                                           Map<String, Map<String, Double>> geneToSampleToValue,
                                           List<String> signature) {
        double total = 0.0d;
        int count = 0;
        for (String gene : signature) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            if (sampleToValue == null)
                continue;
            Double value = sampleToValue.get(sample);
            if (value == null)
                continue;
            count ++;
            total += value;
        }
        return total / count;
    }
    
}
