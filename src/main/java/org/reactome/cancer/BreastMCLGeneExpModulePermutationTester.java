/*
 * Created on Jun 7, 2011
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.EntrezGeneAnalyzer;
import org.reactome.r3.fi.SurvivalAnalysisResult;
import org.reactome.r3.graph.GraphComponent;
import org.reactome.r3.graph.GraphComponentSearchEngine;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.JStatCoxRegressionWrapper;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.ProcessRunner;
import org.reactome.r3.util.R3Constants;


/**
 * This class is used to do some permutation test for breast MCL gene expression modules.
 * @author wgm
 *
 */
public class BreastMCLGeneExpModulePermutationTester {
    private final String DIR_NAME = R3Constants.BREAST_DIR;
//    private final String DIR_NAME = OvarianMCLGeneExpModuleAnalyzer.OVARIAN_DIR_NAME;
    private FileUtility fu = new FileUtility();
    private final String tmp = "tmp";
    private final String tmpFileName = "tmp/TmpModuleToGeneExp.txt";
    
    public BreastMCLGeneExpModulePermutationTester() {
    }
    
    /**
     * THis method is used to check overlapping between greedy components that have better performance than MCL
     * modules and module2.
     * @throws Exception
     */
    @Test
    public void checkOverlapBetweenModule2AndGreedyComponents() throws Exception {
//        String fileName = DIR_NAME + "BetterComponetsFromGeneBasedComponents_JStat_040312.txt";
        String dirName = DIR_NAME + "NEJMGreedySearchComp/";
        String fileName = dirName + "BetterComponetsFromGeneBasedComponents_JStat_040912.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<Integer> selectedIds = new ArrayList<Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            selectedIds.add(new Integer(tokens[0]));
        }
        fu.close();
        
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> module2 = getModule2();
        
        BreastCancerComponentAnalyzer componentHelper = new BreastCancerComponentAnalyzer();
        String compFileName = dirName + "NEJMGeneBasedComponents_JStat_040312.txt";
        List<GraphComponent> components = componentHelper.loadComponentsInOrder(compFileName);
        for (GraphComponent comp : components) {
            if (!selectedIds.contains(comp.getId()))
                continue;
            Set<String> shared = InteractionUtilities.getShared(module2, comp.getAllNodes());
            double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(),
                                                                        comp.getAllNodes().size(),
                                                                        module2.size(),
                                                                        shared.size());
            System.out.println(comp.getId() + "\t" + comp.getAllNodes().size() + "\t" + 
                               shared.size() + "\t" + pvalue);
            System.out.println(shared);
        }
    }
    
    /**
     * This method is used to test the components from greedy search.
     * @throws Exception
     */
    @Test
    public void testGreedyComponents() throws Exception {
        BreastCancerComponentAnalyzer componentHelper = new BreastCancerComponentAnalyzer();
//        String compFileName = DIR_NAME + "GeneBasedComponents_JStat_040312.txt";
//        String dirName = DIR_NAME + "NEJMGreedySearchComp/";
//        String compFileName = dirName + "NEJMGeneBasedComponents_JStat_040312.txt";
        String dirName = DIR_NAME + "GSE4922GreedySearchComp/";
        String compFileName = dirName + "GSE4922GeneBasedComponents_JStat_040312.txt";
        List<GraphComponent> components = componentHelper.loadComponentsInOrder(compFileName);
        for (Iterator<GraphComponent> it = components.iterator(); it.hasNext();) {
            GraphComponent comp = it.next();
            if (comp.getAllNodes().size() < 8 || // Same size as in MCL
                comp.getScore() < 1.30103) {// -Math.log10(0.05)) 
                it.remove(); 
            }
        }
        // Try to filter away redundant components
        componentHelper.filterRedundantComponents(components);
        // Just convert them into gene names
        List<Set<String>> clusters = new ArrayList<Set<String>>();
        int index = 0;
        for (GraphComponent comp : components) {
            clusters.add(comp.getAllNodes());
            System.out.println(index + "\t" + comp.getId() + "\t" + comp.getAllNodes().size());
            index ++;
        }
//        if (true)
//            return;
        String[] srcFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        String[] datasetNames = getDataSetNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        Map<GraphComponent, double[]> compToPvalues = new HashMap<GraphComponent, double[]>();
        for (GraphComponent comp : components) {
            double[] pvalues = new double[srcFileNames.length];
            compToPvalues.put(comp, pvalues);
        }
        for (int i = 0; i < srcFileNames.length; i++) {
            String expFileName = DIR_NAME + srcFileNames[i];
            String clinFileName = DIR_NAME + clinFileNames[i];
//            String tmpFileName = DIR_NAME + datasetNames[i] + "SampleToGreedyComponents_040412.txt";
//            String tmpFileName = dirName + datasetNames[i] + "SampleToFilteredGreedyComponents_040412.txt";
//            String tmpFileName = dirName + datasetNames[i] + "SampleToFilteredGreedyComponents_040512.txt";
            String tmpFileName = dirName + datasetNames[i] + "SampleToFilteredGreedyComponents_040912.txt";
            String[] output = doSurvivalAnalsysi(arrayHelper, 
                                                  clusters, 
                                                  tmpFileName,
                                                  clinFileName, 
                                                  expFileName);
            System.out.println();
            Map<Object, Double> moduleToPvalue = parseSurvivalAnalysisResultForPvalue(output, 
                                                                                      components);
            for (Object obj : moduleToPvalue.keySet()) {
                GraphComponent comp = (GraphComponent) obj;
                double[] pvalues = compToPvalues.get(comp);
                pvalues[i] = moduleToPvalue.get(obj);
            }
        }
        // Search for better results
        double trainingPvalue = 1.75e-10;
        double maxPvalue = 1.64e-4;
        for (Iterator<GraphComponent> it = compToPvalues.keySet().iterator(); it.hasNext();) {
            double[] pvalues = compToPvalues.get(it.next());
            if (pvalues[0] > trainingPvalue)
                it.remove();
            else {
                for (int i = 0; i < pvalues.length; i++) {
                    if (pvalues[i] > maxPvalue) {
                        it.remove();
                        break;
                    }
                }
            }
        }
        System.out.println("\nFiltered components:");
        for (GraphComponent comp : compToPvalues.keySet()) {
            System.out.println(comp.getId() + "\t" + comp.getAllNodes().size() + "\t" + comp.getAllNodes());
        }
    }
    
    /**
     * This method is used to do component search based on Treyer's original algorithm by using
     * CoxPH p-values as scores.
     * @throws Exception
     */
    @Test
    public void doComponentSearch() throws Exception {
        int index = 4;
        String expFileName = getSourceFileNames()[index];
        final String clinFileName = DIR_NAME + getClinFileNames()[index];
        final JStatCoxRegressionWrapper coxregression = new JStatCoxRegressionWrapper();
        final Map<String, Double> sampleToEvent = coxregression.loadSampleToEvent(clinFileName, "OSEVENT");
        final Map<String, Double> sampleToTime = coxregression.loadSampleToSurvival(clinFileName, "OSDURATION");
        
        // Just use the first Nejm data set
        final CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        final Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(DIR_NAME + expFileName);
        GraphComponent.ScoreCalculator scoreCalculate = new GraphComponent.ScoreCalculator() {
            
            @Override
            public double calculateScore(GraphComponent comp) {
//                // Test code below
//                if (true) {
//                    return Math.random();
//                }
                Map<String, Double> sampleToValue = new HashMap<String, Double>();
                Map<String, Integer> sampleToCounts = new HashMap<String, Integer>();
                for (String gene : comp.getAllNodes()) {
                    Map<String, Double> sampleToValue1 = geneToSampleToValue.get(gene);
                    if (sampleToValue1 == null)
                        continue;
                    for (String sample : sampleToValue1.keySet()) {
                        Double value1 = sampleToValue1.get(sample);
                        if (value1 == null)
                            continue;
                        Double value = sampleToValue.get(sample);
                        if (value == null) {
                            sampleToValue.put(sample, value1);
                            sampleToCounts.put(sample, 1);
                        }
                        else {
                            sampleToValue.put(sample, value + value1);
                            sampleToCounts.put(sample, sampleToCounts.get(sample) + 1);
                        }
                    }
                }
                if (sampleToValue.size() > 0) {
                    for (String sample : sampleToValue.keySet()) {
                        Integer counts = sampleToCounts.get(sample);
                        sampleToValue.put(sample, sampleToValue.get(sample) / counts);
                    }
                    double pvalue = coxregression.doSurvivalAnalysis(sampleToValue, sampleToTime, sampleToEvent);
                    if (pvalue == 0.0)
                        return Double.MAX_VALUE;
//                    try {
//                        arrayHelper.generateSampleToGeneExpValue(sampleToValue,
//                                                                 tmpFileName,
//                                                                 "TempCompoennt");
//                        String[] outputs = doSurvivalAnalysis(new File(tmpFileName),
//                                                              new File(clinFileName),
//                                                              "coxph");
//                        String[] tokens = outputs[0].split("\n")[1].split("\t");
//                        Double pvalue = new Double(tokens[2]);
//                        System.out.println("Component " + comp.getAllNodes().size() + ": " + pvalue);
                        return -Math.log10(pvalue);
//                    }
//                    catch(Exception e) {
//                        e.printStackTrace();
//                    }
                }
                return 0;
            }
        };
        long time1 = System.currentTimeMillis();
        GraphComponentSearchEngine searchEngine = new GraphComponentSearchEngine();
        searchEngine.setScoreCalculator(scoreCalculate);
        searchEngine.setMaxDepth(2);
        searchEngine.setSearchDepth(1);
        searchEngine.search(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        List<GraphComponent> components = searchEngine.getFoundComponents();
        System.out.println("Found components: " + components.size());
        BreastCancerComponentAnalyzer componentAnalyzer = new BreastCancerComponentAnalyzer();
//        componentAnalyzer.outputFoundComponents(DIR_NAME + "GeneBasedComponents_JStat_040312.txt",
//                                                components);
//        componentAnalyzer.outputFoundComponents(DIR_NAME + getDataSetNames()[index]+ "GeneBasedComponents_JStat_040312.txt",
//                                                components);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time used: " + (time2 - time1) / 1000 + " seconds");
//        componentAnalyzer.outputFoundComponents(DIR_NAME + "GeneBasedComponents_040312.txt",
//                                                components);
    }
    
    /**
     * This method is used to check overlapping of selected genes from superpc using different
     * data set as the training data set.
     * @throws Exception
     */
    @Test
    public void checkOverlappingFromGeneBasedSuperpc() throws Exception {
        String[] datasetNames = getDataSetNames();
        List<Set<String>> selectedGenes = new ArrayList<Set<String>>();
        Set<String> allGenes = new HashSet<String>();
        for (int i = 0; i < datasetNames.length; i++) {
            System.out.println(datasetNames[i]);
            String subdirName = datasetNames[i] + "ForSingleGenes";
            // Get the list of genes first
//            String posFix = "Top1000Genes_040212.txt";
            String posFix = "SortedGenes_040912.txt";
            fu.setInput(DIR_NAME + subdirName + "/" + datasetNames[i] + posFix);
            String line = fu.readLine();
            String[] genes = line.split("\t");
            for (int j = 1; j < genes.length; j++)
                allGenes.add(genes[j]);
            fu.close();
//            File file = new File(DIR_NAME + subdirName, 
//                                 "SuperpcResults_040212.txt");
            File file = new File(DIR_NAME + subdirName, 
                                 "SuperpcResults_040912.txt");
            fu.setInput(file.getAbsolutePath());
            Set<String> selected = new HashSet<String>();
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("[1] \"Selected gene indices")) {
                    int index = line.indexOf(", ");
                    String sub = line.substring(index + 1, line.length() - 1).trim();
                    String[] tokens = sub.split(", ");
                    for (String token : tokens)
                        selected.add(genes[new Integer(token)]);
                    break;
                }
            }
            fu.close();
            System.out.println("Total genes: " + allGenes.size());
            System.out.println("Selected genes: " + selected.size());
            selectedGenes.add(selected);
        }
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes: " + totalGenes.size());
        System.out.println();
        allGenes = totalGenes; 
        for (int i = 0; i < selectedGenes.size() - 1; i++) {
            Set<String> selected1 = selectedGenes.get(i);
            System.out.println(datasetNames[i]);
            for (int j = i + 1; j < selectedGenes.size(); j ++) {
                Set<String> selected2 = selectedGenes.get(j);
                System.out.println("\t" + datasetNames[j]);
                Set<String> shared = InteractionUtilities.getShared(selected1, selected2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(), 
                                                                            selected1.size(), 
                                                                            selected2.size(), 
                                                                            shared.size());
                System.out.println("\t" + shared.size() + "\t" + pvalue);
            }
        }
        // Check overlapping with module 2 genes
        Set<String> module2Genes = getModule2();
        module2Genes.retainAll(allGenes);
        System.out.println("\nCheck sharing with module2 genes:");
        System.out.println("Module 2 genes: " + module2Genes.size());
        for (int i = 0; i < selectedGenes.size(); i++) {
            Set<String> selected = selectedGenes.get(i);
            Set<String> shared = InteractionUtilities.getShared(selected, module2Genes);
            double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(),
                                                                        module2Genes.size(),
                                                                        selected.size(),
                                                                        shared.size());
            System.out.println(datasetNames[i] + "\t" + shared.size() + "\t" + pvalue);
        }
        // Do a pathway annotations
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        System.out.println("\nPathway annotations");
        for (int i = 0; i < selectedGenes.size(); i++) {
            Set<String> selected = selectedGenes.get(i);
            System.out.println("\n" + datasetNames[i] + " (" + selected.size() + " genes)");
            List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(selected, AnnotationType.Pathway);
            for (GeneSetAnnotation annotation : annotations) {
                if (annotation.getFdr().equals("1.00e+00"))
                    continue;
                if (annotation.getFdr().startsWith("<"))
                    System.out.println(annotation.getTopic() + "\t" + 
                                       annotation.getHitNumber() + "\t" + 
                                       annotation.getFdr());
                else {
                    Double value = new Double(annotation.getFdr());
                    if (value <= 0.05)
                        System.out.println(annotation.getTopic() + "\t" +
                                           annotation.getHitNumber() + "\t" + 
                                           annotation.getFdr());
                }
            }
        }
    }
    
    /**
     * This method is used to generate single gene based files for superpc analysis in order
     * to show our MCL module based method is better than simple gene based method.
     * @throws Exception
     */
    @Test
    public void generateGeneBasedSuperpcFiles() throws Exception {
        String[] srcFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        String[] datasetNames = getDataSetNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        // The following statements are used to pull out genes shared in all data files to avoid any problem in
        // superpc analysis
        Set<String> sharedGenes = new HashSet<String>();
        for (int i = 0; i < datasetNames.length; i++) {
            System.out.println("Source file: " + datasetNames[i]);
            String expFileName = DIR_NAME + srcFileNames[i];
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
            System.out.println("Total genes: " + geneToSampleToValue.size());
            arrayHelper.filterOutGenesWithNullValues(geneToSampleToValue);
            System.out.println("After filtering out genes with null values: " + geneToSampleToValue.size());
            if (sharedGenes.size() == 0)
                sharedGenes.addAll(geneToSampleToValue.keySet());
            else
                sharedGenes.retainAll(geneToSampleToValue.keySet());
            System.out.println("Shared gens: " + sharedGenes.size());
        }
        System.out.println("Total shared genes in all data sets: " + sharedGenes.size());
        List<String> geneList = new ArrayList<String>(sharedGenes);
        Collections.sort(geneList);
        for (int i = 0; i < datasetNames.length; i++) {
            String expFileName = DIR_NAME + srcFileNames[i];
            String clinFileName = DIR_NAME + clinFileNames[i];
            System.out.println("\nSource file: " + expFileName);
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
//            for (int j = 0; j < geneList.size(); j++)
//                System.out.println(j + "\t" + geneList.get(j));
            arrayHelper.generateSampleToGeneExpSubset(geneToSampleToValue,
                                                      tmpFileName,
                                                      geneList);
            File file = new File(tmpFileName);
            String[] output = doSurvivalAnalysis(file,
                                                 new File(clinFileName), 
                                                 "coxph");
//            System.out.println(output[0]);
            final Map<Object, Double> geneToPvalue = parseSurvivalAnalysisResultForPvalue(output,
                                                                                          geneList);
            // Sort genes based on p-values from CoxPH analysis
            Collections.sort(geneList, new Comparator<String>() {
               public int compare(String gene1, String gene2) {
                   Double pvalue1 = geneToPvalue.get(gene1);
                   Double pvalue2 = geneToPvalue.get(gene2);
                   return pvalue1.compareTo(pvalue2);
               }
            });
            String subdirName = datasetNames[i] + "ForSingleGenes";
            File subDir = new File(DIR_NAME, subdirName);
            subDir.mkdir();
            for (int j = 0; j < srcFileNames.length; j++) {
                // Need to reload
                geneToSampleToValue = arrayHelper.loadGeneExp(DIR_NAME + srcFileNames[j]);
//                String tmpFileName = DIR_NAME + subdirName + "/" + datasetNames[j] + "Top1000Genes_040212.txt";
                String tmpFileName = DIR_NAME + subdirName + "/" + datasetNames[j] + "SortedGenes_040912.txt";
                arrayHelper.generateSampleToGeneExpSubset(geneToSampleToValue, 
                                                          tmpFileName, 
                                                          geneList);
            }
        }
    }
    
    @Test
    public void checkOverlappingGenesBasedOnSuperpcFromCrossCheck() throws Exception {
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        String[] geneExpFileNames = getSourceFileNames();
        String[] datasetNames = getDataSetNames();
        List<Set<String>> superpcGenesets = new ArrayList<Set<String>>();
        for (int i = 0; i < datasetNames.length; i++) {
            int index = geneExpFileNames[i].indexOf("/");
//            String clusterFileName = DIR_NAME + "MCL_Clusters_032812_" + geneExpFileNames[i].substring(index + 1);
            String clusterFileName = DIR_NAME + "MCL_Clusters_040612_" + geneExpFileNames[i].substring(index + 1);
            List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
            // Need to get selected clusters from superpc
            String subdirName = datasetNames[i] + "ForMCLTrain";
//            String resultFileName = DIR_NAME + subdirName + "/SuperpcResults_032912.txt";
            String resultFileName = DIR_NAME + subdirName + "/SuperpcResults_040612.txt";
            fu.setInput(resultFileName);
            String line = null;
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("Chosen modules from below")) {
                    index = line.indexOf(":");
                    line = line.substring(index + 1).trim();
                    break;
                }
                else if (line.contains("[1] \"Selected gene indices")) {
                    index = line.indexOf(", ");
                    line = line.substring(index + 1, line.length() - 1).trim();
                    break;
                }
            }
            fu.close();
            String[] tokens = line.split(", ");
            Set<String> genes = new HashSet<String>();
            for (String token : tokens) {
                genes.addAll(clusters.get(new Integer(token) - 1)); // Module indices are from R
            }
            System.out.println("Total genes for " + datasetNames[i] + ": " + genes.size());
            superpcGenesets.add(genes);
        }
        // Check overlapping
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes: " + totalGenes.size());
        System.out.println();
        for (int i = 0; i < datasetNames.length - 1; i++) {
            Set<String> set1 = superpcGenesets.get(i);
            for (int j = i + 1; j < datasetNames.length; j++) {
                Set<String> set2 = superpcGenesets.get(j);
                Set<String> shared = InteractionUtilities.getShared(set1, set2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(),
                                                                            set1.size(),
                                                                            set2.size(),
                                                                            shared.size());
                System.out.println(datasetNames[i] + "\t" + datasetNames[j] + "\t" +
                                   shared.size() + "\t" + pvalue);
            }
        }
    }
    
    /**
     * This method is uesd to check overlapping among significant modules across five breast cancer
     * data sets after using different data set as the training data set.
     * @throws Exception
     */
    @Test
    public void checkOverlappingGenesFromCrossCheck() throws Exception {
        String[] geneExpFileNames = getSourceFileNames();
        // The following is extracted from file MCLCrossCheckResults_032812.xlsx
        int[] significantModules = new int[] {1, 1, 35, 4, 1};
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> sigMclModules = new ArrayList<Set<String>>();
        List<List<Set<String>>> allModules = new ArrayList<List<Set<String>>>();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            int index = geneExpFileNames[i].indexOf("/");
            String clusterFileName = DIR_NAME + "MCL_Clusters_032812_" + geneExpFileNames[i].substring(index + 1);
            List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
            sigMclModules.add(clusters.get(significantModules[i]));
            allModules.add(clusters);
        }
        for (int i = 0; i < sigMclModules.size() - 1; i++) {
            Set<String> module1 = sigMclModules.get(i);
            System.out.println(geneExpFileNames[i]);
            for (int j = i + 1; j < sigMclModules.size(); j++) {
                Set<String> module2 = sigMclModules.get(j);
                Set<String> shared = InteractionUtilities.getShared(module1, module2);
                System.out.println("\t" + geneExpFileNames[j]);
                System.out.println("\tShared: " + shared.size());
            }
        }
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> module2 = getModule2();
        for (int i = 0; i < allModules.size(); i++) {
            List<Set<String>> modules = allModules.get(i);
            System.out.println("\n" + geneExpFileNames[i]);
            for (int j = 0; j < modules.size(); j++) {
                Set<String> module = modules.get(j);
                Set<String> shared = InteractionUtilities.getShared(module2, module);
                double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes.size(),
                                                                            module2.size(),
                                                                            module.size(),
                                                                            shared.size());
                System.out.println(j + "\t" + module.size() + "\t" + 
                                   shared.size() + "\t" + pvalue);
            }
        }
    }
    
    /**
     * In this method, a data set is picked up from 5 breast cancer datasets, and use other four data sets
     * as validation data sets.
     */
    @Test
    public void crossCheckAnalysis() throws Exception {
        long time1 = System.currentTimeMillis();
        String[] geneExpFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        String[] datasetNames = getDataSetNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            // Training data file
            System.out.println("Training file: " + geneExpFileNames[i]);
            String fileName = DIR_NAME + geneExpFileNames[i];
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(fileName);
            Set<String> fisWithCorrs = arrayHelper.calculateGeneExpCorrForFIs(geneToSampleToValue,
                                                                              fis,
                                                                              true,
                                                                              null);
//            // Test code
//            Set<String> fisWithCorrs = new HashSet<String>();
//            for (String fi : fis) {
//                fisWithCorrs.add(fi + "\t1.0");
//            }
            long time3 = System.currentTimeMillis();
            MCLClusterWrapper mclHelper = new MCLClusterWrapper();
            List<Set<String>> clusters = mclHelper.mclCluster(fisWithCorrs, 
                                                              fis, 
                                                              8, 
                                                              0.25d);
            long time4 = System.currentTimeMillis();
            System.out.println("Time for MCL clustering: " + (time4 - time3));
            System.out.println("Clustering results after filtering on size >= 8, correlations >= 0.25:");
            Set<String> totalGenes = new HashSet<String>();
            for (int j = 0; j < clusters.size(); j++) {
                System.out.println(j + "\t" + clusters.get(j).size());
                totalGenes.addAll(clusters.get(j));
            }
            System.out.println("Total genes: " + totalGenes.size());
            if (true)
                return;
            int index = geneExpFileNames[i].lastIndexOf("/");
//            String clusterFileName = DIR_NAME + "MCL_Clusters_032812_" + geneExpFileNames[i].substring(index + 1);
//            String clusterFileName = DIR_NAME + "MCL_Clusters_040412_SingleWeight_" + geneExpFileNames[i].substring(index + 1);
            String clusterFileName = DIR_NAME + "MCL_Clusters_040612_" + geneExpFileNames[i].substring(index + 1);
            clusterHelper.outputNetworkClusters(clusters, clusterFileName);
//            File dir = new File(DIR_NAME, datasetNames[i] + "ForMCLTrain");
            File dir = new File(DIR_NAME, "SingleFIWeight");
            dir.mkdir();
            for (int j = 0; j < geneExpFileNames.length; j ++) {
                String expFileName = DIR_NAME + geneExpFileNames[j];
                String clinFileName = DIR_NAME + clinFileNames[j];
//                File tmpFile = new File(dir, datasetNames[j] + "SampleToMCLModule_032812.txt");
//                File tmpFile = new File(dir, datasetNames[j] + "SampleToMCLModule_040412.txt");
//                File tmpFile = new File(dir, datasetNames[j] + "SampleToMCLModule_040612.txt");
                File tmpFile = new File(dir, datasetNames[j] + "SampleToMCLModule_040912.txt");
                doSurvivalAnalsysi(arrayHelper, 
                                   clusters, 
                                   tmpFile.getAbsolutePath(),
//                                   tmpFileName,
                                   clinFileName, 
                                   expFileName);
            }
            System.out.println();
            break;
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time: " + (time2 - time1));
    }
    
    /**
     * This method is to check the distribution of top genes in module 2 based on survival analysis.
     * @throws Exception
     */
    @Test
    public void checkModule2InGeneBasedSurvivalAnalysis() throws Exception {
//        // The following statements are used to generate a single gene based univariate coxph results
////        String srcFileName = DIR_NAME + "NejmLogRatioNormZScore.txt";
//        String srcFileName = DIR_NAME + "NejmLogRatioNormGlobalZScore_070111.txt";
//        String clinFileName = DIR_NAME + "Nejm_Clin_Simple.txt";
//        // Want to work on this part gene only
//        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
//        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
//        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
//        Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(srcFileName);
//        List<String> fiGeneList = new ArrayList<String>(fiGenes);
//        fiGeneList.retainAll(geneToSampleToValue.keySet());
//        Collections.sort(fiGeneList);
//        for (int i = 0; i < fiGeneList.size(); i++)
//            System.out.println(i + "\t" + fiGeneList.get(i));
//        arrayHelper.generateSampleToGeneExpSubset(geneToSampleToValue,
//                                                  tmpFileName,
//                                                  fiGeneList);
//        File file = new File(tmpFileName);
//        String[] output = doSurvivalAnalysis(file,
//                                             new File(clinFileName), 
//                                             "coxph");
//        System.out.println(output[0]);
        // Use the file compiled from the above result to do enrichment analysis.
        // Genes in this file has been ordered based on p-value
//        String srcFileName = DIR_NAME + "UnivariateCoxphAnalysisOfNejmGenes_060911.txt";
        String srcFileName = DIR_NAME + "UnivariateCoxphAnalysisOfNejmGenes_070511.txt";
        Set<String> module2 = getModule2();
        // Load genes from the source file
        List<String> allGenes = new ArrayList<String>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            allGenes.add(tokens[1]); // The first is gene index based on alphabetical order
        }
        fu.close();
        System.out.println("Total genes: " + allGenes.size());
        // Get the top 100 genes
        List<String> top100Genes = allGenes.subList(0, 100);
        System.out.println("Top 100 genes: " + top100Genes.size());
        // Get shared genes in the top 100 genes
        Set<String> shared = InteractionUtilities.getShared(top100Genes, module2);
        System.out.println("Module 2: " + module2.size());
        System.out.println("Shared: " + shared.size());
        Double pvalue = MathUtilities.calculateHypergeometricPValue(allGenes.size(),
                                                                    top100Genes.size(),
                                                                    module2.size(),
                                                                    shared.size());
        System.out.println("Pvalue from hypergeometric: " + pvalue);
    }
    
    /**
     * This method is used to do survival analysis for a subset genes in a module.
     * @throws Exception
     */
    @Test
    public void doSurvivalAnalysisForSubsetGeneInModule()  throws Exception {
        String[] geneTexts = new String[] {
//                // The following set of genes having p-value < 1.0e-4 from univariate coxph analysis.
//                "BIRC5,AURKB,KIF20A,NCAPD2,CENPN,BUB1,SPC25,RCC2,CDCA8,MAD2L1,NSUN2,SEC13,CCDC99,NDC80,BUB1B",
//                // The following set of genes having p-value < 1.0e-7 from univariate coxph analysis
//                "BIRC5,AURKB,KIF20A,NCAPD2,CENPN,BUB1,SPC25",
//                // The following genes having p-values < 1.0e-9
//                "BIRC5,AURKB,KIF20A",
//                // AURKA has been reported has good prediction power across breast cancer data set.
//                "AURKA",
                // PCNA has been used as another common target for cell proliferation
                "PCNA"
        };
        String[] srcFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        for (String geneText : geneTexts) {
            System.out.println("Genes: " + geneText);
            String[] tokens = geneText.split(",");
            List<Set<String>> clusters = new ArrayList<Set<String>>();
            Set<String> cluster = new HashSet<String>();
            for (String token : tokens) {
                cluster.add(token);
            }
            clusters.add(cluster);
            for (int i = 0; i < srcFileNames.length; i++) {
                String expFileName = DIR_NAME + srcFileNames[i];
                String clinFileName = DIR_NAME + clinFileNames[i];
                System.out.println("Source file: " + expFileName);
                Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
                arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                            tmpFileName,
                                                            clusters);
                File file = new File(tmpFileName);
                String[] output = doSurvivalAnalysis(file,
                                                     new File(clinFileName), 
                                                     "coxph");
                System.out.println(output[0]);
            }
        }
    }
    
    /**
     * This method is used to check the correlation between single gene expression in module
     * 2 and survival times in order to prove that a module is more efficient as a survival 
     * signature.
     * @throws Exception
     */
    @Test
    public void checkSingleGeneInModule2() throws Exception {
        Set<String> module2 = getModule2();
        List<String> module2List = new ArrayList<String>(module2);
        Collections.sort(module2List);
        System.out.println("Second module: " + module2List.size());
        for (int i = 0; i < module2List.size(); i++)
            System.out.println(i + "\t" + module2List.get(i));
        // Gene expression source files
        String[] srcFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        for (int i = 0; i < srcFileNames.length; i++) {
            String expFileName = DIR_NAME + srcFileNames[i];
            String clinFileName = DIR_NAME + clinFileNames[i];
            System.out.println("\nSource file: " + expFileName);
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
            List<String> geneListCopy = new ArrayList<String>(module2List);
            geneListCopy.retainAll(geneToSampleToValue.keySet());
            for (int j = 0; j < geneListCopy.size(); j++)
                System.out.println(j + "\t" + geneListCopy.get(j));
            arrayHelper.generateSampleToGeneExpSubset(geneToSampleToValue,
                                                      tmpFileName,
                                                      geneListCopy);
            File file = new File(tmpFileName);
            String[] output = doSurvivalAnalysis(file,
                                                 new File(clinFileName), 
                                                 "coxph");
            System.out.println(output[0]);
        }
    }

    protected Set<String> getModule2() throws IOException {
        // Get the second module
//        String clusterFileName = R3Constants.BREAST_DIR + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        String clusterFileName = R3Constants.BREAST_DIR + "MCLModulesNejmLogRatioNormGlobalZScore_070111.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
        Set<String> module2 = clusters.get(1);
        return module2;
    }

    protected String[] getClinFileNames() {
        String[] clinFileNames = new String[] {
               "Nejm_Clin_Simple.txt",
               "GSE4922_Clin_Simple.txt",
               "GSE3143/GSE3143_Survival_Info.txt",
//               "GSE1992/GSE1992_OS_Info.txt",
               "GSE18229/GSE18229_Clin_info_Simple_Headers_With_SubTypes.txt",
               "GSE1456/GSE1456_OS_Info.txt"
        };
        return clinFileNames;
    }
    
    protected String[] getDataSetNames() {
        String[] names = new String[] {
                "NEJM",
                "GSE4922",
                "GSE3143",
                "GSE18229",
                "GSE1456"
        };
        return names;
    }

    protected String[] getSourceFileNames() {
        String[] srcFileNames = new String[] {
//                "NejmLogRatioNormZScore.txt",
                "NejmLogRatioNormGlobalZScore_070111.txt",
                "GSE4922FilteredOnSamplesZScore_091409.txt",
//                "GSE3143/GSE3143_MappedGenes_z_on_samples_012111.txt",
                "GSE3143/GSE3143_MappedGenes_z_012111.txt",
//                "GSE1992/GSE1992_Gene_Exp_z_020711.txt",
//                "GSE1992/GSE1992_Gene_Exp_z_080311.txt",
//                "GSE18229/GSE18229_Gene_Exp_z_080311.txt",
                "GSE18229/GSE18229_Gene_Exp_z_111811.txt",
                "GSE1456/GSE1456_z_020911.txt"
        };
        return srcFileNames;
    }
    
    /**
     * Permutating genes in an array data set to generate a random data for MCL clustering.
     * @throws IOException
     */
    @Test
    public void runGeneBasedPermutationTests() throws Exception {
        for (int i = 1; i < 11; i++) {
            String fileName = DIR_NAME + "Nejm/NejmLogRatioNormZScore_random_060711_" + i + ".txt";
            runGeneBasedPermutationTest(fileName);
        }
    }
    
    protected void runGeneBasedPermutationTest(String fileName) throws IOException, Exception {
        System.out.println("Random file: " + fileName + "\n");
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(fileName);
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        // The file size should be around 3.0 M.
        Set<String> fisWithCorrs = arrayHelper.calculateGeneExpCorrForFIs(geneToSampleToValue,
                                                                          fis,
                                                                          true,
                                                                          null);
        MCLClusterWrapper mclHelper = new MCLClusterWrapper();
        List<Set<String>> clusters = mclHelper.mclCluster(fisWithCorrs, 
                                                          fisWithCorrs, 
                                                          8, 
                                                          0.25d);
        System.out.println("Clustering results after filtering on size >= 8, correlations >= 0.25:");
        for (int i = 0; i < clusters.size(); i++) {
            System.out.println(i + "\t" + clusters.get(i).size());
        }
        String tmpFileName = tmp + "/TmpModuleToGeneExp.txt";
        arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                    tmpFileName,
                                                    clusters);
        File file = new File(tmpFileName);
        String clinFileName = "/Users/wgm/datasets/BreastCancer/Nejm_Clin_Simple.txt";
        String[] output = doSurvivalAnalysis(file,
                                             new File(clinFileName), 
                                             "coxph");
        System.out.println("\n" + clinFileName);
        System.out.println(output[0]);
        // The following is based on validation data set
        String[] validationFileNames = new String[] {
                "GSE4922FilteredOnSamplesZScore_091409.txt",
                "GSE3143/GSE3143_MappedGenes_z_on_samples_012111.txt",
                "GSE1992/GSE1992_Gene_Exp_z_020711.txt",
                "GSE1456/GSE1456_z_020911.txt"
        };
        String[] clinFileNames = new String[] {
                "GSE4922_Clin_Simple.txt",
                "GSE3143/GSE3143_Survival_Info.txt",
                "GSE1992/GSE1992_OS_Info.txt",
                "GSE1456/GSE1456_OS_Info.txt"
        };
        for (int i = 0; i < validationFileNames.length; i++) {
            String expFileName = DIR_NAME + validationFileNames[i];
            clinFileName = DIR_NAME + clinFileNames[i];
            doSurvivalAnalsysi(arrayHelper, 
                               clusters, 
                               tmpFileName,
                               clinFileName, 
                               expFileName);
        }
    }

    private String[] doSurvivalAnalsysi(CancerGeneExpressionCommon arrayHelper,
                                    List<Set<String>> clusters,
                                    String tmpFileName, String clinFileName,
                                    String expFileName) throws IOException {
        Map<String, Map<String, Double>> geneToSampleToValue;
        File file;
        String[] output;
        geneToSampleToValue = arrayHelper.loadGeneExp(expFileName);
        System.out.println("Total genes: " + geneToSampleToValue.size());
        arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                    tmpFileName,
                                                    clusters);
        file = new File(tmpFileName);
        output = doSurvivalAnalysis(file,
                                    new File(clinFileName), 
                                    "coxph");
        System.out.println("\n" + clinFileName);
        System.out.println(output[0]);
        return output;
    }

    
    private String[] doSurvivalAnalysis(File scoreFile,
                                        File clinFile,
                                        String model) throws IOException {
        return doSurvivalAnalysis(scoreFile, 
                                  clinFile,
                                  model,
                                  null);
    }
    
    protected String[] doSurvivalAnalysis(File scoreFile,
                                        File clinFile,
                                        String model,
                                        Integer moduleIndex) throws IOException {
        ProcessRunner runner = new ProcessRunner();
        List<String> parameters = new ArrayList<String>();
        // TODO: Have to make sure this is in the path. Otherwise it should be changed
        parameters.add("Rscript");
        parameters.add(R3Constants.survivalScript);
        parameters.add(scoreFile.getAbsolutePath());
        parameters.add(clinFile.getAbsolutePath());
        parameters.add(model);
        if (moduleIndex != null)
            parameters.add(moduleIndex.toString());
        String[] output = runner.runScript(parameters.toArray(new String[0]));
        SurvivalAnalysisResult result = new SurvivalAnalysisResult();
        result.setOutput(output[0]);
        result.setError(output[1]);
        return output;
    }
    
    @Test
    public void generateRandomizedGeneExpFile() throws IOException {
        String srcFileName = DIR_NAME + "NejmLogRatioNormZScore.txt";
        List<String> genes = new ArrayList<String>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            genes.add(line.substring(0, index));
        }
        fu.close();
        // perform permutations
        Collections.shuffle(genes);
        String targetFileName = DIR_NAME + "Nejm/NejmLogRatioNormZScore_random_060711_10.txt";
        fu.setInput(srcFileName);
        fu.setOutput(targetFileName);
        line = fu.readLine();
        fu.printLine(line);
        int geneIndex = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String gene = genes.get(geneIndex);
            geneIndex ++;
            fu.printLine(gene + line.substring(index));
        }
        fu.close();
    }
    
    /**
     * Use this method to check the performance of other signatures in these five breast cancer data sets.
     * @throws Exception
     */
    @Test
    public void checkOtherSignature() throws Exception {
//        Set<String> signatures = new BreastMCLGeneExpModuleAnalyzer().loadOncotypeDXGenes(false);
//        Set<String> signatures = fu.loadInteractions(DIR_NAME + "nejm_table2/Gene70_Partial.txt");
//        Set<String> signatures = new BreastMCLGeneExpModuleAnalyzer().loadMCLModule2();
        Set<String> signatures = fu.loadInteractions(TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt");
        System.out.println("Total signature genes: " + signatures.size());
        List<String> list1 = new ArrayList<String>(signatures);
        Collections.sort(list1);
        System.out.println(list1);
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        boolean needEdit = signatures.contains("AKAP2");
        signatures = entrezAnalyzer.normalizeGeneNames(signatures);
        // A manual fix for NEJM gene signature: some mess in the original NCBI gene name file
        if (needEdit && signatures.contains("PALM2-AKAP2")) {
            signatures.remove("PALM2-AKAP2");
            signatures.add("AKAP2");
        }
        System.out.println("After name normalization: " + signatures.size());
        list1 = new ArrayList<String>(signatures);
        Collections.sort(list1);
        System.out.println(list1);
        
        List<Set<String>> list = new ArrayList<Set<String>>();
        list.add(signatures);
//        String[] geneExpFileNames = new String[] {
//                "NejmLogRatioNormGlobalZScore_070111.txt",
//                "GSE4922FilteredOnSamplesZScore_091409.txt",
//                "GSE3143/GSE3143_MappedGenes_z_012111.txt",
//                "GSE18229/GSE18229_Gene_Exp_z_080311.txt",
////                "GSE18229/GSE18229_Gene_Exp_z_111811.txt",
//                "GSE1456/GSE1456_z_020911.txt"
//        };
//        String[] clinFileNames = new String[] {
//                "Nejm_Clin_Simple.txt", // For the training data set
//                "GSE4922_Clin_Simple.txt",
//                "GSE3143/GSE3143_Survival_Info.txt",
//                "GSE18229/GSE18229_Clin_info_Simple_Headers.txt",
//                "GSE1456/GSE1456_OS_Info.txt"
//        };
//        String[] geneExpFileNames = new String[] {
//                "TothillDataset/GSE9891_Gene_Exp_111811.txt",
////                "data_072610/TCGA_batch9-15_17-19_21-22_24.UE.No_Normal.txt",
//                "data_published_in_nature/TCGA_batch9-15_17-19_21-22_24.UE.No_Normal.No_Mutated_in_Module6_7.txt",
//                "GSE13876/GSE13876_Exp_log2_norm_gene_avg.txt",
//                "GSE26712/GSE26712_Gene_Exp_111911.txt"
//        };
//        String[] clinFileNames = new String[] {
//                "TothillClinInfo111811.txt",
////                "data_031910/Batches9-22_tcga_OV_clinical_csv.2010-01-25-noDates.txt",
//                "data_published_in_nature/2010-09-11380C-Table_S1.2_with_edit.txt",
//                "GSE13876/GSE13876_Clin.txt",
//                "GSE26712/GSE26712ClinInfo.txt"
//        };
        String[] geneExpFileNames = getSourceFileNames();
        String[] clinFileNames = getClinFileNames();
        List<Map<String, Map<String, Double>>> geneToSampleToValues = new ArrayList<Map<String,Map<String,Double>>>();
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        for (int i = 0; i < geneExpFileNames.length; i++) {
            Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(DIR_NAME + geneExpFileNames[i]);
            geneToSampleToValues.add(geneToSampleToValue);
        }
        for (int i = 0; i < clinFileNames.length; i++) {
            Map<String, Map<String, Double>> geneToSampleToValue = geneToSampleToValues.get(i);
            String tmpFileName = tmp + "/TmpModuleToGeneExp.txt";
            arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                        tmpFileName,
                                                        list);
            File file = new File(tmpFileName);
            String clinFileName = DIR_NAME + clinFileNames[i];
            String[] output = doSurvivalAnalysis(file,
                                                 new File(clinFileName), 
                                                "coxph",
                                                0);
            System.out.println(geneExpFileNames[i]);
            System.out.println(output[0]);
            System.out.println(output[1]);
        }
    }
    
    /**
     * Compare performance with randomly selected gene set without considering the FI network.
     * @throws Exception
     */
    @Test
    public void permutationTestModuleSignficanceUsingGenesOnly() throws Exception {
        long time1 = System.currentTimeMillis();
        // Load the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> allGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Load gene expression data file
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        // The following files are used for test
        String[] validationFileNames = getSourceFileNames();

        RandomData randomizer = new RandomDataImpl();
        // Need to load all expression together to increase the speed
        List<Map<String, Map<String, Double>>> geneToSampleToValueList = new ArrayList<Map<String,Map<String,Double>>>();
        for (int i = 0; i < validationFileNames.length; i++) {
            Map<String, Map<String, Double>> tmp = arrayHelper.loadGeneExp(DIR_NAME + validationFileNames[i]);
            geneToSampleToValueList.add(tmp);
        }
        String[] clinFileNames = getClinFileNames();
        int significantModuleCount = 0;
        List<String> sigModulesInOneRun = new ArrayList<String>();
        int permutationNumber = 1000; 
//        double pvalueCutoff = 0.05;
        int geneSize = getModule2().size();
        double pvalueCutoff = 1.64E-04;
        System.out.println("Target gene size: " + geneSize);
        System.out.println("Target p-value: " + pvalueCutoff);
        List<Set<String>> modules = new ArrayList<Set<String>>();
        for (int i = 0; i < permutationNumber; i++) {
            System.out.println("Permutation " + i + "...");
            Set<String> randomGenes = MathUtilities.randomSampling(allGenes, geneSize, randomizer);
            modules.clear();
            modules.add(randomGenes);
            sigModulesInOneRun.clear();
            for (int j = 0; j < clinFileNames.length; j++) {
                Map<String, Map<String, Double>> geneToSampleToValue = geneToSampleToValueList.get(j);
                String tmpFileName = tmp + "/TmpModuleToGeneExppermutationTestModuleSignficanceUsingGenesOnly.txt";
                arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                            tmpFileName,
                                                            modules);
                File file = new File(tmpFileName);
                String clinFileName = DIR_NAME + clinFileNames[j];
                String[] output = doSurvivalAnalysis(file,
                                                     new File(clinFileName), 
                                                    "coxph");
//                System.out.println(output[0]);
                List<String> sigModules = parseSurvivalAnalysisResult(output,
                                                                      pvalueCutoff);
                if (j == 0) {
                    sigModulesInOneRun.addAll(sigModules);
                }
                else
                    sigModulesInOneRun.retainAll(sigModules);
                if (sigModulesInOneRun.size() == 0)
                    break;
            }
            System.out.println("Significant modules: " + sigModulesInOneRun.size());
            significantModuleCount += sigModulesInOneRun.size();
        }
        System.out.println("Total significant modules in permutation " + permutationNumber + ": " + significantModuleCount);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
    }
    
    /**
     * This permutation test is used to check if the results from the Nejm data set is significant across
     * several different data sets based on random permutating genes.
     * @throws Exception
     */
    @Test
    public void permutationTestModuleSignificance() throws Exception {
        long time1 = System.currentTimeMillis();
        // Load the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        // Load gene expression data file
        CancerGeneExpressionCommon arrayHelper = new CancerGeneExpressionCommon();
        // The following files are used for test
        String[] validationFileNames = getSourceFileNames();
        String srcFileName = DIR_NAME + validationFileNames[0];
        Map<String, Map<String, Double>> geneToSampleToValue = arrayHelper.loadGeneExp(srcFileName);
//        Map<String, Map<String, Double>> originalGeneToSampleValue = new HashMap<String, Map<String, Double>>(geneToSampleToValue);
        // Create two lists for permutation
        List<String> geneList = new ArrayList<String>(geneToSampleToValue.keySet());
        List<Map<String, Double>> sampleToValueList = new ArrayList<Map<String, Double>>(geneToSampleToValue.values());
        
        Map<String, Map<String, Double>> geneToSampleToValueRandom = new HashMap<String, Map<String,Double>>();
        RandomData randomizer = new RandomDataImpl();

        // Need to load all expression together to increase the speed
        List<Map<String, Map<String, Double>>> validationGeneToSampleToValue = new ArrayList<Map<String,Map<String,Double>>>();
        for (int i = 1; i < validationFileNames.length; i++) {
            Map<String, Map<String, Double>> tmp = arrayHelper.loadGeneExp(DIR_NAME + validationFileNames[i]);
            validationGeneToSampleToValue.add(tmp);
        }
        String[] clinFileNames = getClinFileNames();
//        List<Map<String, Double>> sampleToEvents = new ArrayList<Map<String, Double>>();
//        List<Map<String, Double>> sampleToTimes = new ArrayList<Map<String, Double>>();
//        JStatCoxRegressionWrapper coxWrapper = new JStatCoxRegressionWrapper();
//        for (int i = 0; i < clinFileNames.length; i++) {
//            Map<String, Double> sampleToEvent = coxWrapper.loadSampleToEvent(DIR_NAME + clinFileNames[i], 
//                                                                             "OSEVENT");
//            sampleToEvents.add(sampleToEvent);
//            Map<String, Double> sampleToTime = coxWrapper.loadSampleToSurvival(DIR_NAME + clinFileNames[i],
//                                                                               "OSDURATION");
//            sampleToTimes.add(sampleToTime);
//        }
        int significantModuleCount = 0;
        List<String> sigModulesInOneRun = new ArrayList<String>();
        int permutationNumber = 1000; 
//        double pvalueCutoff = 0.05;
//        double pvalueCutoff = 1.99E-04;
        double pvalueCutoff = 1.64E-04; // As of April 6, 2012
        System.out.println("Pvalue cutoff: " + pvalueCutoff);
        MCLClusterWrapper mclCluster = new MCLClusterWrapper();
        for (int i = 0; i < permutationNumber; i++) {
            System.out.println("Permutation " + i + "...");
            Object[] geneListRandom = randomizer.nextSample(geneList, geneList.size());
            Object[] sampleToValueRandom = randomizer.nextSample(sampleToValueList, sampleToValueList.size());
            geneToSampleToValueRandom.clear();
            for (int j = 0; j < geneListRandom.length; j++) {
                String gene = (String) geneListRandom[j];
                Map<String, Double> sampleToValue = (Map<String, Double>) sampleToValueRandom[j];
                geneToSampleToValueRandom.put(gene, sampleToValue);
            }
//            // The following is used for test the code
//            geneToSampleToValueRandom = originalGeneToSampleValue;
//            Map<String, Double> sampleToValue = geneToSampleToValueRandom.get("FLJ20011");
//            String sample = sampleToValue.keySet().iterator().next();
//            Double value = sampleToValue.get(sample);
//            System.out.println(sample + ": " + value);
            
            Set<String> fisWithCorrs = arrayHelper.calculateGeneExpCorrForFIs(geneToSampleToValueRandom,
                                                                              fis,
                                                                              true,
                                                                              null);
            List<Set<String>> modules = mclCluster.mclCluster(fisWithCorrs, 
                                                              fis,
                                                              8, 
                                                              0.25d);
            // Search for module that is significant across all samples
            sigModulesInOneRun.clear();
            for (int j = 0; j < clinFileNames.length; j++) {
                if (j == 0)
                    geneToSampleToValue = geneToSampleToValueRandom;
                else
                    geneToSampleToValue = validationGeneToSampleToValue.get(j - 1);
                String tmpFileName = tmp + "/TmpModuleToGeneExpForpermutationTestModuleSignificance.txt";
                arrayHelper.generateSampleToGeneExpClusters(geneToSampleToValue,
                                                            tmpFileName,
                                                            modules);
                File file = new File(tmpFileName);
                String clinFileName = DIR_NAME + clinFileNames[j];
                String[] output = doSurvivalAnalysis(file,
                                                     new File(clinFileName), 
                                                    "coxph");
                List<String> sigModules = parseSurvivalAnalysisResult(output,
                                                                      pvalueCutoff);
                // Test code to use JStatSoft package. However, it is easy to get a Matrix is singular error.
//                List<Map<String, Double>> sampleToValues = arrayHelper.generateSampleToGeneExpForClusters(geneToSampleToValue,
//                                                                                                          modules);
//                List<Double> pvalues = coxWrapper.doSurvivalAnalyses(sampleToValues, 
//                                                                     sampleToTimes.get(j),
//                                                                     sampleToEvents.get(j));
//                List<String> sigModules = new ArrayList<String>();
//                for (int k = 0; k < pvalues.size(); k++) {
//                    if (pvalues.get(k) <= pvalueCutoff) {
//                        sigModules.add(k + 1 + "");
//                    }
//                }
                if (j == 0) {
//                    System.out.println(output[0]);
                    sigModulesInOneRun.addAll(sigModules);
                }
                else
                    sigModulesInOneRun.retainAll(sigModules);
                if (sigModulesInOneRun.size() == 0)
                    break;
            }
            System.out.println("Significant modules: " + sigModulesInOneRun.size());
            significantModuleCount += sigModulesInOneRun.size();
        }
        System.out.println("Total significant modules in permutation " + permutationNumber + ": " + significantModuleCount);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
    }
    
    /**
     * A helper method to get module index which is significant
     * @param output
     * @return
     * @throws ParseException 
     */
    private List<String> parseSurvivalAnalysisResult(String[] output,
                                                     double cutoff) throws ParseException  {
        List<String> rtn = new ArrayList<String>();
        // The first element is the result
        String result = output[0];
        String[] lines = result.split("\n");
        // The first line is the header
        for (int i = 1; i < lines.length; i++) {
            String[] tokens = lines[i].split("\t");
            double pvalue = new Double(tokens[2]);
            if (pvalue <= cutoff)
                rtn.add(tokens[0]);
        }
        return rtn;
    }
    
    private Map<Object, Double> parseSurvivalAnalysisResultForPvalue(String[] output,
                                                                      List<?> geneList) throws ParseException {
        Map<Object, Double> geneToPvalue = new HashMap<Object, Double>();
        String[] lines = output[0].split("\n");
        for (int i = 1; i < lines.length; i++) {
            String[] tokens = lines[i].split("\t");
            Object gene = geneList.get(new Integer(tokens[0]));
            geneToPvalue.put(gene,
                              new Double(tokens[2]));
        }
        return geneToPvalue;
    }
    
    @Test
    public void testParseNumber() throws ParseException {
        String text = "1.751864e-10";
        Double value = new Double(text);
        double dValue = value.doubleValue();
        System.out.println("Value: " + value);
        System.out.println("In number: " + dValue);
//        value = 1.99E-04d;
    }
    
}
