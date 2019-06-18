/*
 * Created on May 18, 2006
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.Value;

// Cannot find jsmile for Bayesian Network in maven. Disable the code for the time being
//import smile.Network;
//import smile.Network.BayesianAlgorithmType;
//import smile.learning.DataSet;
//import smile.learning.EM;
//import smile.learning.NaiveBayes;
//import smile.learning.TextParser;


/**
 * This class is used for analysis related to Bayesian Network.
 * @author guanming
 *
 */
public class BayesianNetworkAnalyzer extends TestCase {
    // This map is used to hold all values to be used for training
    private Map<String, Value> values;
    // Cut-off value
    private double CUT_OFF = 0.51;
    // Use to control if GO data should be discretized.
    private boolean DISCRETIZE_GO = false;
    
    public BayesianNetworkAnalyzer() {
        values = new HashMap<String, Value>();
    }
    
//    /**
//     * Use thie method to create training dataset for SimpleBayesianNetwork.
//     * @throws Exception
//     */
//    public void generateTrainingDataSet() throws Exception {
//        String posDataFileName = "results/interaction/ReactomeInteractions.txt";
//        String negDataFileName = "results/interaction/NoInteractionsForTrain.txt";
//        prepareDataSet(posDataFileName, negDataFileName, true);
//        // Localization data cannot be used in the training since negative training dataset
//        // is created based on localization. This way makes all colocalization protein pairs
//        // are functional interactions.
//        //generateLocalization(posIntSet, negIntSet);
//        exportData("results/TrainingData.txt");
//    }
//    
//    public void generateTrainingDataSetForNonY2HInteractions() throws Exception {
//        String posDataFileName = "results/interaction/HumanNonY2H.txt";
//        String negDataFileName = "results/interaction/HumanNonY2HNoInteractions.txt";
//        prepareDataSet(posDataFileName, negDataFileName, false);
//        exportDataForWEKA("results/TrainingDataForNonY2H.arff");
//    }
//    
//    public void generateTestDataSet() throws Exception {
//        String posDataFileName = "results/interaction/PantherInteractions.txt";
//        String negDataFileName = "results/interaction/NoInteractionsForTest.txt";
//        prepareDataSet(posDataFileName, negDataFileName, true);
//        exportData("results/TestingData.txt");
//    }
//    
//    /**
//     * Use this method to create training dataset for Naive Bayesian Classifier.
//     * @throws Exception
//     */
//    public void generateTrainingDataSetForNBC() throws Exception {
//        String posDataFileName = "results/interaction/ReactomeInteractions.txt";
//        String negDataFileName = "results/interaction/NoInteractionsForTrain.txt";
//        prepareDataSet(posDataFileName, negDataFileName, true);
//        exportDataForNBC();
//    }
//
//    private void prepareDataSet(String posFileName,
//                                String negFileName,
//                                boolean needHumanPPI) throws IOException, Exception {
//        values.clear();
//        Set<String> posIntSet = new HashSet<String>();
//        loadProteinPair(posFileName, posIntSet);
//        // Initialize Values
//        for (String pair : posIntSet) {
//            Value value = getValueObject(pair);
//            value.functionalInteraction = Boolean.TRUE;
//        }
//        Set<String> negIntSet = new HashSet<String>();
//        loadProteinPair(negFileName, negIntSet);
//        for (String pair : negIntSet) {
//            Value value = getValueObject(pair);
//            value.functionalInteraction = Boolean.FALSE;
//        }
//        Set<String> total = new HashSet<String>(posIntSet);
//        total.addAll(negIntSet);
//        System.out.println("total dataset: " + total.size());
//        System.out.println("total positive set: " + posIntSet.size());
//        System.out.println("total negative set: " + negIntSet.size());
//        if (needHumanPPI)
//            generateInteraction(posIntSet, negIntSet);
//        generateOrthoInteraction(posIntSet, negIntSet);
//        //generateGeneExp(posIntSet, negIntSet);
//        //generateGeneExpForIndividualDataset(posIntSet, negIntSet);
//        generateGOData(posIntSet, negIntSet);
//    }
//    
//    private void exportDataForNBC() throws IOException {
//        String fileName = "results/TrainingDataForNBC.txt";
//        StringBuilder builder = new StringBuilder();
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        fu.printLine("FunctionalInteraction\tY2H\tNonY2H" +
//                     "\tGeneExp\tGOBPOccurence\tGOBPDepth" +
//                     "\tGOMFOccurence\tGOMFDepth");
//        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
//            String pair = it.next();
//            Value value = values.get(pair);
//            builder.append(value.functionalInteraction).append("\t");
//            if (value.y2h == null)
//                builder.append("false\t");
//            else
//                builder.append(value.y2h).append("\t");
//            if (value.nonY2h == null)
//                builder.append("false\t");
//            else
//                builder.append(value.nonY2h).append("\t");
//            builder.append(value.geneExp).append("\t");
//            if (value.goBPOccurence == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goBPOccurence).append("\t"); 
//            if (value.goBPDepth == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goBPDepth).append("\t");
//            if (value.goMFOccurence == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goMFOccurence).append("\t");
//            if (value.goMFDepth == null)
//                builder.append("*");
//            else
//                builder.append("s").append(value.goMFDepth);
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
//    }
//    
//    private void exportData(String fileName) throws IOException {
//        StringBuilder builder = new StringBuilder();
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        fu.printLine("FunctionalInteraction\tPhysicalInteraction\tY2H\tNonY2H" +
//                     "\tGeneExp\tGOBP\tGOBPOccurence\tGOBPDepth\tGOMF" +
//                     "\tGOMFOccurence\tGOMFDepth");
//        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
//            String pair = it.next();
//            Value value = values.get(pair);
//            builder.append(value.functionalInteraction).append("\t");
//            builder.append("*\t"); // unknown
//            if (value.y2h == null)
//                builder.append("false\t");
//            else
//                builder.append(value.y2h).append("\t");
//            if (value.nonY2h == null)
//                builder.append("false\t");
//            else
//                builder.append(value.nonY2h).append("\t");
//            builder.append(value.geneExp).append("\t");
//            // for GOBP
//            builder.append("*\t");
//            if (value.goBPOccurence == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goBPOccurence).append("\t"); 
//            if (value.goBPDepth == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goBPDepth).append("\t");
//            // for GOMF
//            builder.append("*\t");
//            if (value.goMFOccurence == null)
//                builder.append("*\t");
//            else
//                builder.append("s").append(value.goMFOccurence).append("\t");
//            if (value.goMFDepth == null)
//                builder.append("*");
//            else
//                builder.append("s").append(value.goMFDepth);
//            // localization
//            //if (value.localization == null)
//            //    builder.append("*");
//            //else
//            //    builder.append(value.localization);
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
//    }
//    
////    private void generateLocalization(Set<String> positiveSet,
////                                      Set<String> negativeSet) throws IOException {
////        // Load protein to localization map
////        GODataAnalyzer analyzer = new GODataAnalyzer();
////        Map<String, List<String>> protein2Localization = analyzer.loadProtein2Locations();
////        int c = 0;
////        int index = 0;
////        Value valueObj;
////        for (String pair : positiveSet) {
////            valueObj = getValueObject(pair);
////            index = pair.indexOf(" ");
////            String id1 = pair.substring(0, index);
////            String id2 = pair.substring(index + 1);
////            List<String> loc1 = protein2Localization.get(id1);
////            List<String> loc2 = protein2Localization.get(id2);
////            // These should be treated as "unknown"
////            if (loc1 == null || loc2 == null)
////                continue;
////            if (analyzer.isColocalized(loc1, loc2)) {
////                c ++;
////                valueObj.localization = Boolean.TRUE;
////            }
////            else
////                valueObj.localization = Boolean.FALSE;
////        }
////        System.out.println("Positive Col: " + c);
////        c = 0;
////        for (String pair : negativeSet) {
////            valueObj = getValueObject(pair);
////            index = pair.indexOf(" ");
////            String id1 = pair.substring(0, index);
////            String id2 = pair.substring(index + 1);
////            List<String> loc1 = protein2Localization.get(id1);
////            List<String> loc2 = protein2Localization.get(id2);
////            if (analyzer.isColocalized(loc1, loc2)) {
////                c ++;
////                valueObj.localization = Boolean.TRUE;
////            }
////            else
////                valueObj.localization = Boolean.FALSE;
////        }
////        System.out.println("Negative Col: " + c);
////    }
//    
//    private void generateGOData(Set<String> posIntSet,
//                                Set<String> negIntSet) throws IOException {
//        // Need to load annoations
//        GODataAnalyzer analyzer = new GODataAnalyzer();
//        Map<String, Set<String>> prot2BPMap = new HashMap<String, Set<String>>();
//        Map<String, Set<String>> prot2MFMap = new HashMap<String, Set<String>>();
//        analyzer.generateProtein2GOMap(prot2MFMap, prot2BPMap);
//        // First: BP Occurence
//        generateGOOccurenceData(true,
//                                prot2BPMap,
//                                posIntSet,
//                                negIntSet,
//                                analyzer);
//        // Second: MF Occurence
//        generateGOOccurenceData(false,
//                                prot2MFMap,
//                                posIntSet,
//                                negIntSet,
//                                analyzer);
//        // Third: BP Depth
//        Map<String, Integer> depthMap = new HashMap<String, Integer>();
//        // Load depth values
//        FileUtility fu = new FileUtility();
//        fu.setInput("results/go/GODepth.txt");
//        String line = null;
//        int index = 0;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            depthMap.put(line.substring(0, index), Integer.parseInt(line.substring(index + 1)));
//        }
//        fu.close();
//        generateGODepthData(depthMap, prot2BPMap, posIntSet, negIntSet, analyzer, true);
//        // Fourth: MF Depth
//        generateGODepthData(depthMap, prot2MFMap, posIntSet, negIntSet, analyzer, false);
//    }
//    
//    private void generateGODepthData(Map<String, Integer> depthMap,
//                                     Map<String, Set<String>> prot2GOMap,
//                                     Set<String> positiveSet,
//                                     Set<String> negativeSet,
//                                     GODataAnalyzer analyzer,
//                                     boolean isForBP) throws IOException {
//        int index = 0;
//        Value valueObj;
//        int c = 0;
//        for (String pair : positiveSet) {
//            valueObj = getValueObject(pair);
//            index = pair.indexOf(" ");
//            String id1 = pair.substring(0, index);
//            String id2 = pair.substring(index + 1);
//            int depth = analyzer.findSharedDepth(prot2GOMap.get(id1),
//                                                 prot2GOMap.get(id2),
//                                                 depthMap);
//            if (depth > Integer.MIN_VALUE) {
//                if (isForBP)
//                    valueObj.goBPDepth = depth;
//                else
//                    valueObj.goMFDepth = depth;
//                c ++;
//            }
//        }
//        String message = "Shared for ";
//        if (isForBP)
//            message += "GOBPDepth:";
//        else
//            message += "GOMFDepth:";
//        System.out.println("Positive " + message + c);
//        c = 0;
//        for (String pair : negativeSet) {
//            valueObj = getValueObject(pair);
//            index = pair.indexOf(" ");
//            String id1 = pair.substring(0, index);
//            String id2 = pair.substring(index + 1);
//            int depth = analyzer.findSharedDepth(prot2GOMap.get(id1),
//                                                 prot2GOMap.get(id2),
//                                                 depthMap);
//            if (depth > Integer.MIN_VALUE) {
//                if (isForBP)
//                    valueObj.goBPDepth = depth;
//                else
//                    valueObj.goMFDepth = depth;
//                c ++;
//            }
//        }
//        System.out.println("Negative " + message + c);
//    }
//    
//    private List<int[]> loadOccurenceCategories(String fileName) throws IOException {
//        List<int[]> categories = new ArrayList<int[]>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        int index = 0;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf(" ");
//            int[] bounds = new int[] {
//                    Integer.parseInt(line.substring(0, index)),
//                    Integer.parseInt(line.substring(index + 1))
//            };
//            categories.add(bounds);
//        }
//        return categories;
//    }
//    
//    private void generateGOOccurenceData(boolean isForBP,
//                                         Map<String, Set<String>> prot2GOMap,
//                                         Set<String> positiveSet, 
//                                         Set<String> negativeSet,
//                                         GODataAnalyzer analyzer) throws IOException {
//        String fileName, categoryFileName;
//        if (isForBP) {
//            fileName = "results/go/GOPOccurence.txt";
//            categoryFileName = "results/go/GOBPOccurenceCategories.txt";
//        }
//        else {
//            fileName = "results/go/GOFOccurence.txt";
//            categoryFileName = "results/go/GOMFOccurenceCategories.txt";
//        }
//        // Load categories
//        List<int[]> categories = loadOccurenceCategories(categoryFileName);
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        Map<String, Integer> occurenceMap = new HashMap<String, Integer>();
//        String line = null;
//        int index = 0;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            occurenceMap.put(line.substring(0, index), 
//                             Integer.parseInt(line.substring(index + 1)));
//        }
//        fu.close();
//        int c = 0;
//        Value valueObj;
//        for (String pair : positiveSet) {
//            valueObj = getValueObject(pair);
//            index = pair.indexOf(" ");
//            String id1 = pair.substring(0, index);
//            String id2 = pair.substring(index + 1);
//            int occurence = analyzer.findSharedTermOccurence(prot2GOMap.get(id1),
//                                                             prot2GOMap.get(id2),
//                                                             occurenceMap);
//            if (occurence < Integer.MAX_VALUE) {
//                if (DISCRETIZE_GO) {
//                    int category = getOccurenceCategory(occurence, categories);
//                    if (isForBP)
//                        valueObj.goBPOccurence = category;
//                    else
//                        valueObj.goMFOccurence = category;
//                }
//                else {
//                    if (isForBP)
//                        valueObj.goBPOccurence = occurence;
//                    else
//                        valueObj.goMFOccurence = occurence;
//                }
//                c++;
//            }
//        }
//        String message = "Shared for ";
//        if (isForBP)
//            message += "GOBPOccurence:";
//        else
//            message += "GOMFOccurence:";
//        System.out.println("Positive " + message + c);
//        c = 0;
//        for (String pair : negativeSet) {
//            valueObj = getValueObject(pair);
//            index = pair.indexOf(" ");
//            String id1 = pair.substring(0, index);
//            String id2 = pair.substring(index + 1);
//            int occurence = analyzer.findSharedTermOccurence(prot2GOMap.get(id1),
//                                                             prot2GOMap.get(id2),
//                                                             occurenceMap);
//            if (occurence < Integer.MAX_VALUE) {  
//                if (DISCRETIZE_GO) {
//                    int category = getOccurenceCategory(occurence, categories);
//                    if (isForBP)
//                        valueObj.goBPOccurence = category;
//                    else
//                        valueObj.goMFOccurence = category;
//                }
//                else {
//                    if (isForBP)
//                        valueObj.goBPOccurence = occurence;
//                    else
//                        valueObj.goMFOccurence = occurence;
//                }
//                c++;
//            }
//        }
//        System.out.println("Negative " + message + c);
//    }
//    
//    private int getOccurenceCategory(int occurence, List<int[]> categories) {
//        for (int i = 0; i < categories.size(); i++) {
//            int[] bounds = categories.get(i);
//            if (occurence >= bounds[0] && occurence <= bounds[1])
//                return i;
//        }
//        throw new IllegalStateException("Cannot find a category for occurence: " + occurence);
//    }
//    
//    private void generateGeneExp(Set<String> posIntSet, 
//                                 Set<String> negIntSet) throws Exception {
//        // For checking gene expression data
//        Class.forName("com.mysql.jdbc.Driver");
//        Connection connection = DriverManager.getConnection("jdbc:mysql:/localhost:3306/BNDataSource?" +
//                "user=root&password=macmysql01");
//        String query = "SELECT value FROM GeneExpPair WHERE pair=?";
//        PreparedStatement stat = connection.prepareStatement(query);  
//        int posC = 0;
//        int negC = 0;
//        Value valueObj;
//        long time1 = System.currentTimeMillis();
//        //O14569 O60304
//        //String pair = "O14569" + " " + "O60304";
//        //stat.setString(1, pair);
//        for (String pair : posIntSet) { 
//            valueObj = getValueObject(pair);
//            valueObj.geneExp = "non"; // Default
//            stat.setString(1, pair);
//            ResultSet result = stat.executeQuery();
//            if (result.next()) {
//                String value = result.getString(1);
//                if (value.equals("+")) {
//                    posC ++;
//                    valueObj.geneExp = "pos";
//                }
//                else {
//                    negC ++;
//                    valueObj.geneExp = "neg";
//                }
//            }
//            result.close();
//        }
//        System.out.printf("Positive from GeneExp: pos %d, neg %d%n", posC, negC);
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for getting GeneExp for Positive: " + (time2 - time1));
//        posC = negC = 0;
//        for (String pair : negIntSet) { 
//            valueObj = getValueObject(pair);
//            valueObj.geneExp = "non";
//            stat.setString(1, pair);
//            ResultSet result = stat.executeQuery();
//            if (result.next()) {
//                String value = result.getString(1);
//                if (value.equals("+")) {
//                    posC ++;
//                    valueObj.geneExp = "pos";
//                }
//                else {
//                    negC ++;
//                    valueObj.geneExp = "neg";
//                }
//            }
//            result.close();
//        }
//        stat.close();
//        connection.close();
//        System.out.printf("Negative from GeneExp: pos %d, neg %d%n", posC, negC);
//        long time3 = System.currentTimeMillis();
//        System.out.println("Time for getting GeneExp for Negative: " + (time3 - time2));
//    }
//    
//    private void generateGeneExpForIndividualDataset(Set<String> posIntSet,
//                                                     Set<String> negIntSet) throws Exception {
//        // Load MicroArrayAnalyzer
//        MicroarrayDataAnalyzer analyzer = new MicroarrayDataAnalyzer();
//        List<String> gdsList = analyzer.getGDSList();
//        for (String gds : gdsList) {
//            analyzer.loadData("results/microarray/" + gds + ".ave.txt");
//            Value valueObj;
//            long time1 = System.currentTimeMillis();
//            //String pair = "O14569 O60304"; separted by a single space
//            int index = 0;
//            String id1, id2;
//            Float cor = null;
//            for (String pair : posIntSet) { 
//                valueObj = getValueObject(pair);
//                index = pair.indexOf(" ");
//                id1 = pair.substring(0, index);
//                id2 = pair.substring(index + 1);
//                cor = analyzer.calculateCorrelation(id1, id2);
//                valueObj.addGeneExpValue(cor);
//            }
//            long time2 = System.currentTimeMillis();
//            System.out.println("Time for getting GeneExp for Positive: " + (time2 - time1));
//            for (String pair : negIntSet) { 
//                valueObj = getValueObject(pair);
//                index = pair.indexOf(" ");
//                id1 = pair.substring(0, index);
//                id2 = pair.substring(index + 1);
//                cor = analyzer.calculateCorrelation(id1, id2);
//                valueObj.addGeneExpValue(cor);
//            }
//            long time3 = System.currentTimeMillis();
//            System.out.println("Time for getting GeneExp for Negative: " + (time3 - time2));
//            System.out.println("Done..." + gds);
//        }
//    }
//    
//    private void generateOrthoInteraction(Set<String> posIntSet,
//                                         Set<String> negIntSet) throws Exception {
//        Set<String> flyNonY2H = loadOrthoInteraction("results/interaction/flyNonY2HInteractionHsaUni.txt");
//        Set<String> flyY2H = loadOrthoInteraction("results/interaction/flyY2HInteractionHsaUni.txt");
//        Set<String> wormY2H = loadOrthoInteraction("results/interaction/wormY2HInteractionHsaUni.txt");
//        // Check occurences
//        int cFlyNonY2H = 0;
//        int cFlyY2H = 0;
//        int cWormY2H = 0;
//        for (String pair : posIntSet) {
//            Value value = getValueObject(pair);
//            if (flyNonY2H.contains(pair)) {
//                cFlyNonY2H ++;
//                value.flyNonY2H = Boolean.TRUE;
//            }
//            if (flyY2H.contains(pair)) {
//                cFlyY2H ++;
//                value.flyY2H = Boolean.TRUE;
//            }
//            if (wormY2H.contains(pair)) {
//                cWormY2H ++;
//                value.wormY2H = Boolean.TRUE;
//            }
//        }
//        System.out.printf("Positive: FlyNonY2H %d, FlyY2H %d, WormY2H %d%n", cFlyNonY2H, cFlyY2H, cWormY2H);
//        cFlyNonY2H = 0;
//        cFlyY2H = 0;
//        cWormY2H = 0;
//        for (String pair : negIntSet) {
//            Value value = getValueObject(pair);
//            if (flyNonY2H.contains(pair)) {
//                cFlyNonY2H ++;
//                value.flyNonY2H = Boolean.TRUE;
//            }
//            if (flyY2H.contains(pair)) {
//                cFlyY2H ++;
//                value.flyY2H = Boolean.TRUE;
//            }
//            if (wormY2H.contains(pair)) {
//                cWormY2H ++;
//                value.wormY2H = Boolean.TRUE;
//            }
//        }
//        System.out.printf("Negative: FlyNonY2H %d, FlyY2H %d, WormY2H %d%n", cFlyNonY2H, cFlyY2H, cWormY2H);
//    }
//    
//    private Set<String> loadOrthoInteraction(String fileName) throws IOException {
//        Set<String> set = new HashSet<String>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        while ((line = fu.readLine()) != null)
//            set.add(line);
//        fu.close();
//        return set;
//    }
//    
//    private void generateInteraction(Set<String> posIntSet,
//                                     Set<String> negIntSet) throws Exception {
//        // Load Protein Interaction Data
//        Set<String> y2hInteractions = new HashSet<String>();
//        Set<String> nonY2HInteractions = new HashSet<String>();
//        loadInteractions(y2hInteractions, nonY2HInteractions);
//        // Check occurences
//        int cy2=0, cny2=0;
//        for (String pair : posIntSet) {
//            Value value = getValueObject(pair);
//            if (y2hInteractions.contains(pair)) {
//                cy2 ++;
//                value.y2h = Boolean.TRUE;
//            }
//            if (nonY2HInteractions.contains(pair)) {
//                cny2 ++;
//                value.nonY2h = Boolean.TRUE;
//            }
//        }
//        System.out.printf("Positive: Y2H %d, NonY2H %d%n", cy2, cny2);
//        cy2 = 0;
//        cny2 = 0;
//        for (String pair : negIntSet) {
//            Value value = getValueObject(pair);
//            if (y2hInteractions.contains(pair)) {
//                cy2 ++;
//                value.y2h = Boolean.TRUE;
//            }
//            if (nonY2HInteractions.contains(pair)) {
//                cny2 ++;
//                value.nonY2h = Boolean.TRUE;
//            }
//        }
//        System.out.printf("Negative: Y2H %d, NonY2H %d%n", cy2, cny2);
//    }
//    
//    private Value getValueObject(String pair) {
//        Value value = values.get(pair);
//        if (value == null) {
//            value = new Value();
//            values.put(pair, value);
//        }
//        return value;
//    }
//    
//    private void loadProteinPair(String fileName, Set<String> set) throws IOException {
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        while ((line = fu.readLine()) != null) {
//            set.add(line);
//        }
//        fu.close();
//    }
//    
//    private void loadInteractions(Set<String> y2hInteractions, Set<String> nonY2HInteractions)
//                    throws Exception {
//        String interactionFileName = "results/interaction/interactions.xml";
//        SAXBuilder builder = new SAXBuilder();
//        Document document = builder.build(interactionFileName);
//        Element root = document.getRootElement();
//        List intElms = root.getChildren("interaction");
//        System.out.println("total interactions: " + intElms.size());
//        for (Iterator it = intElms.iterator(); it.hasNext();) {
//            Element elm = (Element) it.next();
//            Element expElm = elm.getChild("experimentTypes");
//            String[] expNames = expElm.getText().split(",");
//            Element actorsElm = elm.getChild("interactors");
//            List actorsList = actorsElm.getChildren("interactor");
//            Element tmp = (Element) actorsList.get(0);
//            String id1 = tmp.getAttributeValue("id");
//            tmp = (Element) actorsList.get(1);
//            String id2 = tmp.getAttributeValue("id");
//            String prtPair = createProteinPair(id1, id2);
//            for (String expName : expNames) {
//                if (expName.startsWith("two hybrid"))
//                    y2hInteractions.add(prtPair);
//                else
//                    nonY2HInteractions.add(prtPair);
//            }
//        }
//        System.out.println("Total Y2H: " + y2hInteractions.size());
//        System.out.println("Total Non Y2H: " + nonY2HInteractions.size());
//    }
//    
//    private String createProteinPair(String id1, String id2) {
//        if (id1.compareTo(id2) < 0) 
//            return id1 + " " + id2;
//        return id2 + " " + id1;
//    }
//    
//    public void learnBNParameters() throws Exception {
//        // Set up BN
//        Network network = new Network();
//        network.readFile("results/SimpleFunctionalInteractionBN.xdsl");
//        // Load training data
//        DataSet dataset = parseDataFile("results/TrainingData.txt", network);
//        // Need to map nodes in dataset to nodes in network before learning
//        String[] nodes = network.getAllNodeIds();
//        int size = dataset.getVariableCount();
//        for (int i = 0; i < size; i++) {
//            String name = dataset.getVariableId(i);
//            int handle = network.getNode(name);
//            dataset.setNodeHandle(i, network.getNode(name));
//        }
//        // Use EM to learn parameters
//        long time1 = System.currentTimeMillis();
//        EM em = new EM();
//        em.setEqSampleSize(0);
//        em.setRandomizeParameters(true);
//        em.learn(dataset, network);
//        network.writeFile("results/LearnedSimpleBN.xdsl");
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for learning: " + (time2 - time1));
//    }
//    
//    private DataSet parseDataFile(String fileName, Network network) throws IOException {
//        DataSet dataset = new DataSet();
//        // Need the state names for each node
//        String[] variableNames = network.getAllNodeIds();
//        Map<String, List<String>> node2States = new HashMap<String, List<String>>();
//        for (String name : variableNames) {
//            String[] states = network.getOutcomeIds(name);
//            node2States.put(name, Arrays.asList(states));
//        }
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        // The first line is header
//        String line = fu.readLine();
//        String[] ids = line.split("\t");
//        for (String id : ids)
//            dataset.addIntVariable(id);
//        // After the header, all are data lines
//        int variableCount = ids.length;
//        int index = 0;
//        int record = 0;
//        int stateIndex = 0;
//        while ((line = fu.readLine()) != null) {
//            String[] values = line.split("\t");
//            dataset.addEmptyRecord();
//            for (String v : values) {
//                if (v.equals("*"))
//                    stateIndex = DataSet.DefaultMissingInt;
//                else {
//                    // Need to find out the index of v in the output array
//                    List<String> outputs = node2States.get(ids[index]);
//                    stateIndex = outputs.indexOf(v);
//                }
//                dataset.setInt(index, record, stateIndex);
//                index ++;
//            }
//            record ++;
//            index = 0;
//        }
//        return dataset;
//    }
//    
//    public void testBNInference() throws Exception {
//        Network network = new Network();
//        //network.readFile("results/LearnedSimpleBN.xdsl");
//        network.readFile("results/LearnedNaiveBayes.xdsl");
//        network.setBayesianAlgorithm(BayesianAlgorithmType.Pearl);
//        // update to calculate all values
//        network.updateBeliefs();
//        // Set true for localization as evidence
//        //network.setEvidence("Y2H", "true");
//        //network.setEvidence("NonY2H", "false");
//        network.updateBeliefs();
//        // Get functionalInteraction value
//        double[] values = network.getNodeValue("FunctionalInteraction");
//        for (double v : values)
//            System.out.println("values: " + v);
//    }
//    
//    public void learnNaiveBayesClassifier() throws Exception {
//        TextParser parser = new TextParser();
//        parser.setUseHeader(true);
//        if (parser.parse("results/TrainingDataForNBC.txt") != Network.OK) {
//            System.err.println("Cannot parse the datafile!");
//            System.exit(1);
//        }
//        DataSet ds = parser.getDataSet();
//        NaiveBayes nb = new NaiveBayes();
//        nb.setClassVariableId("FunctionalInteraction");
//        System.out.println("Starting learning...");
//        long time1 = System.currentTimeMillis();
//        Network network = nb.learn(ds);
//        network.writeFile("results/LearnedNaiveBayes.xdsl");
//        System.out.println("Learning is done.");
//        long time2 = System.currentTimeMillis();
//        System.out.println("Time for learning: " + (time2 - time1));
//    }
//    
//    public void testBayesianNetwork() throws Exception {
//        // Load the network
//        Network network = new Network();
//        network.readFile("results/LearnedSimpleBN.xdsl");
//        testBN(network);
//    }
//    
//    private List<Value> loadData(String fileName) throws IOException {
//        List<Value> values = new ArrayList<Value>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = fu.readLine(); // The header
//        while ((line = fu.readLine()) != null) {
//            Value value = new Value();
//            //FunctionalInteraction PhysicalInteraction Y2H NonY2H GeneExp GOBP    
//            //GOBPOccurence GOBPDepth GOMF GOMFOccurence GOMFDepth
//            String[] tokens = line.split("\t");
//            value.functionalInteraction = Boolean.valueOf(tokens[0]);
//            if (!tokens[2].equals("*"))
//                value.y2h = Boolean.valueOf(tokens[2]);
//            if (!tokens[3].equals("*"))
//                value.nonY2h = Boolean.valueOf(tokens[3]);
//            if (!tokens[4].equals("*"))
//                value.geneExp = tokens[4];
//            if (!tokens[6].equals("*"))
//                value.goBPOccurence = Integer.valueOf(tokens[6].substring(1)); // Strip "s"
//            if (!tokens[7].equals("*"))
//                value.goBPDepth = Integer.valueOf(tokens[7].substring(1));
//            if (!tokens[9].equals("*"))
//                value.goMFOccurence = Integer.valueOf(tokens[9].substring(1));
//            if (!tokens[10].equals("*"))
//                value.goMFDepth = Integer.valueOf(tokens[10].substring(1));
//            values.add(value);
//        }
//        fu.close();
//        return values;
//    }
//    
//    private void testBN(Network network) throws Exception {
//        network.setTarget("FunctionalInteraction", true);
//        List<Value> values = loadData("results/TestingData.txt");
//        int posCorrectCount = 0;
//        for (Value value : values) {
//            if (value.functionalInteraction && testInstance(network, value))
//                posCorrectCount ++;
//        }
//        System.out.println("True Positive: " + posCorrectCount);
//        int negCorrectCount = 0;
//        for (Value value : values) {
//            if (!value.functionalInteraction && testInstance(network, value)) {
//                negCorrectCount ++;
//            }
//        }
//        System.out.println("True Negative: " + negCorrectCount);
//    }
//    
//    private boolean testInstance(Network network, Value value) {
//        network.clearAllEvidence();
//        network.updateBeliefs();
//        if (value.y2h != null) {
//            network.setEvidence("Y2H", value.y2h.toString());
//        }
//        if (value.nonY2h != null)
//            network.setEvidence("NonY2H", value.nonY2h.toString());
//        network.setEvidence("GeneExp", value.geneExp.toString());
//        if (value.goBPDepth != null)
//            network.setEvidence("GOBPOccurence", "s" + value.goBPDepth);
//        if (value.goBPOccurence != null)
//            network.setEvidence("GOBPOccurence", "s" + value.goBPOccurence);
//        if (value.goMFDepth != null)
//            network.setEvidence("GOMFDepth", "s" + value.goMFDepth);
//        if (value.goMFOccurence != null)
//            network.setEvidence("GOMFOccurence", "s" + value.goMFOccurence);
//        network.updateBeliefs();
//        // Get the values for FunctionalInteractions
//        double[] values = network.getNodeValue("FunctionalInteraction");
//        String[] outcomes = network.getOutcomeIds("FunctionalInteraction");
//        for (int i = 0; i < outcomes.length; i++) {
//            if (outcomes[i].equals(value.functionalInteraction.toString())) {
//                //System.out.println(outcomes[i] + " " + values[i]);
//                if (values[i] > CUT_OFF)
//                    return true;
//                else
//                    return false;
//            }
//        }
//        return false;
//    }
//    
//    public void testNaiveBayesClassifier() throws Exception {
//        Network network = new Network();
//        network.readFile("results/LearnedNaiveBayes.xdsl");
//        testBN(network);
//    }
//    
//    public void generateTrainingDataForWEKA() throws Exception {
//        String posDataFileName = "results/interaction/ReactomeInteractions.txt";
//        String negDataFileName = "results/interaction/NoInteractionsForTrain.txt";
//        prepareDataSet(posDataFileName, negDataFileName, true);
//        exportDataForWEKA("results/TrainingDataWithGrid.arff");
//    }
//    
//    public void generateTestingDataForWEKA() throws Exception {
//        String posDataFileName = "results/interaction/PantherInteractions.txt";
//        String negDataFileName = "results/interaction/NoInteractionsForTest.txt";
//        prepareDataSet(posDataFileName, negDataFileName, true);
//        exportDataForWEKA("results/TestingDataWithOrtho.arff");
//    }
//    
//    private void exportDataForWEKA(String fileName) throws IOException {
//        StringBuilder builder = new StringBuilder();
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        fu.printLine("@relation FunctionTraining");
//        fu.printLine("");
//        fu.printLine("@attribute FunctionalInteraction {true, false}");
//        fu.printLine("@attribute Y2H {true, false}");
//        fu.printLine("@attribute NonY2H {true, false}");
//        fu.printLine("@attribute FlyY2H {true, false}");
//        fu.printLine("@attribute FlyNonY2H {true, false}");
//        fu.printLine("@attribute WormY2H {true, false}");
//        //MicroarrayDataAnalyzer analyzer = new MicroarrayDataAnalyzer();
//        //List<String> gdsList = analyzer.getGDSList();
//        //for (String gds : gdsList) {
//        //    fu.printLine("@attribute " + gds + " real");
//        //}
//        //fu.printLine("@attribute GeneExp real");
//        fu.printLine("@attribute GeneExp {pos, neg, non}");
//        fu.printLine("@attribute GOBPOccurence real");
//        fu.printLine("@attribute GOBPDepth real");
//        fu.printLine("@attribute GOMFOccurence real");
//        fu.printLine("@attribute GOMFDepth real");
//        fu.printLine("");
//        fu.printLine("@data");
//        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
//            String pair = it.next();
//            Value value = values.get(pair);
//            builder.append(value.functionalInteraction).append(",");
//            if (value.y2h == null)
//                builder.append("false,");
//            else
//                builder.append(value.y2h).append(",");
//            if (value.nonY2h == null)
//                builder.append("false,");
//            else
//                builder.append(value.nonY2h).append(",");
//            if (value.flyY2H == null)
//                builder.append("false,");
//            else
//                builder.append(value.flyY2H).append(",");
//            if (value.flyNonY2H == null)
//                builder.append("false,");
//            else
//                builder.append(value.flyNonY2H).append(",");
//            if (value.wormY2H == null)
//                builder.append("false,");
//            else
//                builder.append(value.wormY2H).append(",");
////            List<Float> geneExpValues = value.geneExpList;
////            for (Float f : geneExpValues) {
////                if (f == null)
////                    builder.append("?,");
////                else
////                    builder.append(f).append(",");
////            }
//            if (value.geneExp == null)
//                builder.append("?,");
//            else
//                builder.append(value.geneExp).append(",");
//            if (value.goBPOccurence == null)
//                builder.append("?,");
//            else
//                builder.append(value.goBPOccurence).append(","); 
//            if (value.goBPDepth == null)
//                builder.append("?,");
//            else
//                builder.append(value.goBPDepth).append(",");
//            if (value.goMFOccurence == null)
//                builder.append("?,");
//            else
//                builder.append(value.goMFOccurence).append(",");
//            if (value.goMFDepth == null)
//                builder.append("?");
//            else
//                builder.append(value.goMFDepth);
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
//    }
//    
//    public void checkoverlapWithPhysicalInteraction() throws Exception {
//        // Load Protein Interaction Data
//        Set<String> y2hInteractions = new HashSet<String>();
//        Set<String> nonY2HInteractions = new HashSet<String>();
//        loadInteractions(y2hInteractions, nonY2HInteractions);
//        Set<String> totalInteractions = new HashSet<String>();
//        totalInteractions.addAll(y2hInteractions);
//        totalInteractions.addAll(nonY2HInteractions);
//        String fileName = "results/microarray/GeneExpFromPavlidis.txt"; 
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        int index = 0;
//        String pair = null;
//        int c = 0;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            pair = line.substring(index + 1);
//            if (totalInteractions.contains(pair))
//                c ++;
//        }
//        System.out.println("Overlapping between interactions and gene expression: " + c);
//    }
//    
//    public void checkOverlapWithReactome() throws Exception {
//        String inFile = "results/interaction/ReactomeInteractions.txt";
//        FileUtility fu = new FileUtility();
//        Set<String> reactomeInteractions = fu.loadInteractions(inFile);
//        Set<String> reactomeIds = InteractionUtilities.grepIDsFromInteractions(reactomeInteractions);
//        inFile = "results/interaction/PantherInteractions.txt";
//        Set<String> pantherInteractions = fu.loadInteractions(inFile);
//        Set<String> pantherIds = InteractionUtilities.grepIDsFromInteractions(pantherInteractions);
//        System.out.println("Total IDs in Panther interactions: " + pantherIds.size());
//        System.out.println("Total Interaction in Panther: " + pantherInteractions.size());
//        pantherInteractions.removeAll(reactomeInteractions);
//        System.out.println("Interactions in Panther but not in Reactome: " + pantherInteractions.size());
//        //pantherIds.addAll(reactomeIds);
//        //System.out.println("Total IDs from both Reactome and Panther interactions: " + pantherIds.size());
//        // Load Protein Interaction Data
//        inFile = "results/interaction/yeastY2HInteractionHsaUni.txt";
//        Set<String> yeastY2H = fu.loadInteractions(inFile);
//        inFile = "results/interaction/yeastNonY2HInteractionHsaUni.txt";
//        Set<String> yeastNonY2H = fu.loadInteractions(inFile);
//        // Check overlapping with Reactome interactions
//        int c = 0;
//        int c1 = 0;
//        for (String pair : reactomeInteractions) {
//            if (yeastY2H.contains(pair))
//                c ++;
//            if (yeastNonY2H.contains(pair))
//                c1 ++;
//        }
//        System.out.printf("Reactome in Yeast Y2H: %s, NonY2H: %s%n", c, c1);
//        c = c1 = 0;
//        for (String pair : pantherInteractions) {
//            if (yeastY2H.contains(pair))
//                c ++;
//            if (yeastNonY2H.contains(pair))
//                c1 ++;
//        }
//        System.out.printf("Panther in Yeast Y2H: %s, NonY2H: %s%n", c, c1);
//    }
//    
//    public void checkHumanInteraction() throws Exception {
//        Set<String> y2hSet = new HashSet<String>();
//        Set<String> nonY2HSet = new HashSet<String>();
//        loadInteractions(y2hSet, nonY2HSet);
//        Set<String> negSet = new HashSet<String>();
//        generateOrthoInteraction(y2hSet, negSet);
//        generateGeneExp(y2hSet, negSet);
//        generateGOData(y2hSet, negSet);
//    }
    
}
