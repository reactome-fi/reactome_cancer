/*
 * Created on May 4, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do SSP based breast cancer data analysis.
 * @author wgm
 *
 */
public class BreastSubtypeSSPAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public BreastSubtypeSSPAnalyzer() {
    }
    
    /**
     * Changed from old name check306GenesSignature to checkOverlapBetweenModule2AndSignature.
     * @throws Exception
     */
    @Test
    public void checkOverlapBetweenModule2AndSignature() throws Exception {
//        Set<String> signatureGenes = loadSSPGenes();
        Set<String> signatureGenes = loadPam50Genes();
        System.out.println("Total genes from 306 signature: " + signatureGenes.size());
        // Load module 2
        // Doing overlapping analysis
        String clusterFileName = R3Constants.BREAST_DIR + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
        //        for (int i = 0; i < clusters.size(); i++)
        //            System.out.println(i + ": " + clusters.get(i).size());
        // Get cluster 6, which is the one we want
        Set<String> module2 = clusters.get(6);
        System.out.println("Second module (" + module2.size() + "): " + module2);
        int totalGene = R3Constants.TOTAL_HUMAN_GENES;
        System.out.println("\nOverlapping analysis:");
        Set<String> shared = InteractionUtilities.getShared(signatureGenes,
                                                            module2);
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGene, 
                                                                    module2.size(),
                                                                    signatureGenes.size(),
                                                                    shared.size());
        System.out.println(shared.size() + "\t" + pvalue);
        System.out.println("Shared genes: " + shared);
    }
    
    /**
     * This method is used to load 50 genes used in PAM50 classifier.
     * @return
     * @throws IOException
     */
    private Set<String> loadPam50Genes() throws IOException {
        String fileName = R3Constants.BREAST_DIR + "PAM50/bioclassifier_R/pam50_centroids.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> rtn = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            rtn.add(tokens[0]);
        }
        return rtn;
    }

    /**
     * This method is used to load SSP genes.
     * @return
     * @throws IOException
     */
    protected Set<String> loadSSPGenes() throws IOException {
        String fileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/306genes-X-249samples-X-5subtypes.txt";
        Map<String, String> ugidToSymbol = loadUGIDToSymbol();
        // Load all genes
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> signatureGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String symbol = ugidToSymbol.get(tokens[0]);
            if (symbol != null)
                signatureGenes.add(symbol);
        }
        fu.close();
        return signatureGenes;
    }
    
    @Test
    public void checkGeneSignatures() throws IOException, MathException {
        String fileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/Gene_Signatures_SSP_UNC_Mapped.txt";
        Map<String, List<String>> clusterToGenes = new HashMap<String, List<String>>();
        Map<String, List<String>> clusterToSymbols = new HashMap<String, List<String>>();
        fu.setInput(fileName);
        String line = null;
        List<String> genes = null;
        List<String> symbols = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#")) {
                // Get the cluster name
                String clusterName = line.substring(1).trim();
                genes = new ArrayList<String>();
                symbols = new ArrayList<String>();
                clusterToGenes.put(clusterName, genes);
                clusterToSymbols.put(clusterName, symbols);
            }
            else {
                String[] tokens = line.split("\t");
                genes.add(tokens[0]);
                if (tokens.length > 1 && tokens[1].length() > 0)
                    symbols.add(tokens[1]);
            }
        }
        fu.close();
        // Check the total gene names
        int total = 0;
        int symbolTotal = 0;
        for (String cluster : clusterToGenes.keySet()) {
            genes = clusterToGenes.get(cluster);
            symbols = clusterToSymbols.get(cluster);
            System.out.println(cluster + ": " + genes.size());
            System.out.println(cluster + " symbols: " + symbols.size());
            total += genes.size();
            symbolTotal += symbols.size();
        }
        System.out.println("Total genes: " + total);
        System.out.println("Total symbols: " + symbolTotal);
        // Doing overlapping analysis
        String clusterFileName = R3Constants.BREAST_DIR + "MCL_Clusters_I_50_FIsWithNejmGeneExpAbsCorr_011111.txt";
        NetworkClusterAnalyzer clusterHelper = new NetworkClusterAnalyzer();
        List<Set<String>> clusters = clusterHelper.loadNetworkClusters(clusterFileName);
//        for (int i = 0; i < clusters.size(); i++)
//            System.out.println(i + ": " + clusters.get(i).size());
        // Get cluster 6, which is the one we want
        Set<String> module2 = clusters.get(6);
        System.out.println("Second module (" + module2.size() + "): " + module2);
        int totalGene = R3Constants.TOTAL_HUMAN_GENES;
        System.out.println("\nOverlapping analysis:");
        Set<String> shared = null;
        for (String cluster : clusterToSymbols.keySet()) {
            symbols = clusterToSymbols.get(cluster);
            shared = InteractionUtilities.getShared(symbols, module2);
            double pvalue = MathUtilities.calculateHypergeometricPValue(totalGene, 
                                                                        module2.size(),
                                                                        symbols.size(),
                                                                        shared.size());
            System.out.println(cluster + "\t" + symbols.size() + "\t" + shared.size() + "\t" + pvalue);
        }
        // Get the shared for the last genes
        System.out.println("Shared genes: " + InteractionUtilities.joinStringElements(", ", shared));
//        Map<String, String> geneNameToSymbol = loadGeneNameToSymbol();
//        List<String> symbols = new ArrayList<String>();
//        System.out.println("Mapped to symbols...");
//        // The following code is used to generate a mapping from gene names to symbols
//        for (String cluster : clusterToGenes.keySet()) {
//            genes = clusterToGenes.get(cluster);
//            Collections.sort(genes);
//            // Check how many can be mapped
//            System.out.println("# " + cluster);
//            for (String name : genes) {
//                String symbol = searchSymbolOnName(name, geneNameToSymbol);
//                System.out.println(name + "\t" + (symbol == null ? "" : symbol));
//                if (symbol == null) {
////                    System.err.println(name + " cannot be mapped!");
//                }
//                else {
//                    symbols.add(symbol);
//                    //System.out.println(name + "\t" + symbol);
//                }
//            }
////            System.out.println(cluster + ": " + symbols.size());
//        }
//        System.out.println("Total mapped: " + symbols.size());
    }
    
    private String searchSymbolOnName(String name, Map<String, String> nameToSymbol) {
        String symbol = nameToSymbol.get(name);
        if (symbol != null)
            return symbol;
        // Search based on the first match
        for (String tmpName : nameToSymbol.keySet()) {
            symbol = nameToSymbol.get(tmpName);
            tmpName = tmpName.replaceAll("\\(|\\)", "");
            if (tmpName.equals(name))
                return symbol;
            if (tmpName.startsWith(name))
                return symbol;
        }
        return null;
    }
    
    Map<String, String> loadGeneNameToSymbol() throws IOException {
        String fileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/Hu-et-al-Intrinsic-List-with-Annotation.txt";
        Map<String, String> nameToSymbol = new HashMap<String, String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            String[] tokens = line.split("\t");
            if (tokens.length < 4)
                continue; // Data not found
            String name = tokens[2].trim();
            // In case there is any quotation marks
            name = name.replace("\"", "");
            String symbol = tokens[3].trim();
            if (symbol.length() == 0)
                continue;
            // Check if a mapping has been there
            String oldSymbol = nameToSymbol.get(name);
            if (oldSymbol == null)
                nameToSymbol.put(name, symbol);
            else if (!oldSymbol.equals(symbol)) {
                throw new IllegalStateException(name + " has more than one gene symbol!");
            }
        }
        fu.close();
        return nameToSymbol;
    }
    
    private Map<String, String> loadUGIDToSymbol() throws IOException {
        Map<String, String> ugidToSymbol = new HashMap<String, String>();
        String fileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/Hu-et-al-Intrinsic-List-with-Annotation.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            String[] tokens = line.split("\t");
            if (tokens.length < 4)
                continue; // Data not found
            String symbol = tokens[3].trim();
            if (symbol.length() == 0)
                continue;
            ugidToSymbol.put(tokens[0].trim(), symbol);
        }
        fu.close();
        return ugidToSymbol;
    }
    
    private Map<String, String> loadGeneNameToUGID() throws IOException {
        String mapFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE1456GenesToCLID.txt";
        Map<String, String> nameToUGID1 = generateNameToUGID(mapFileName);
        mapFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE4922GenesToCLID.txt";
        Map<String, String> nameToUGID2 = generateNameToUGID(mapFileName);
        System.out.println("Size of map1: " + nameToUGID1.size());
        System.out.println("Size of map2: " + nameToUGID2.size());
        Map<String, String> nameToUGID = new HashMap<String, String>(nameToUGID1);
        nameToUGID.putAll(nameToUGID2);
        System.out.println("Merged: " + nameToUGID.size());
        return nameToUGID;
    }
    
    
    
    @Test
    public void mergeSubTypesToNejmClinFile() throws Exception {
//        String srcFileName = R3Constants.BREAST_DIR + "Nejm_ClinicalData_Table.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "Nejm_Full_Clinical_Table_With_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/NejmSubTypesOnSSPNoDWD_NormZScore.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "Nejm_Full_Clinical_Table_With_SubTypes.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "Nejm/Nejm_Full_Clinical_Table_With_Pam50_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "Nejm/pam50_classifier/Nejm_Sub_Type_Pam50.txt";
        
        // For GSE4922
//        String srcFileName = R3Constants.BREAST_DIR + "GSE4922_Clin_Full.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE4922_Clin_Full_With_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE4922SubTypesOnSSPNoDWD.txt";
      String srcFileName = R3Constants.BREAST_DIR + "GSE4922_Clin_Full_With_SubTypes.txt";
      String targetFileName = R3Constants.BREAST_DIR + "GSE4922/GSE4922_Clin_Full_With_PAM50_SubTypes.txt";
      String subTypeFile = R3Constants.BREAST_DIR + "GSE4922/pam50_classifier/GSE4922_PAM50_SubTypes.txt";


        // For GSE1456
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_OS_Info.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_OS_Info_With_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE1456SubTypesOnSSPNoDWD.txt";
//      String srcFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_OS_Info_With_SubTypes.txt";
//      String targetFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_OS_Info_With_Pam50_SubTypes.txt";
//      String subTypeFile = R3Constants.BREAST_DIR + "GSE1456/pam50_classifier/GSE1456_PAM50_SubTypes.txt";

        
        // GSE3143
//        String srcFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_Survival_Info.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_Survival_Info_With_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE3143SubTypesOnSSPNoDWD.txt";
//      String srcFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_Survival_Info_With_SubTypes.txt";
//      String targetFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_Survival_Info_With_Pam50_SubTypes.txt";
//      String subTypeFile = R3Constants.BREAST_DIR + "GSE3143/pam50_classifier/GSE3143_PAM50_SubTypes.txt";

        
        // GSE1992
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1992/GSE1992_Clin_Info_Simple_Header.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE1992/GSE1992_Clin_Info_Simple_Header_With_SubTypes.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE1992SubTypesOnSSPNoDWD.txt";
        
////        GSE18229
//        String srcFileName = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Clin_info_Simple_Headers.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Clin_info_Simple_Headers_With_SubTypes.txt";
////        String subTypeFile = R3Constants.BREAST_DIR + "GSE18229/GSE18229SubTypesOnSSPNoDWD.txt";
//        String subTypeFile = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Sub_Type_Pam50_Self.txt";
        
        Map<String, String> sampleToType = new HashMap<String, String>();
        fu.setInput(subTypeFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            sampleToType.put(tokens[0], tokens[2]);
        }
        fu.close();
        fu.setInput(srcFileName);
        fu.setOutput(targetFileName);
        line = fu.readLine();
        fu.printLine(line + "\tSub_Type");
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            // Only NEJM samples need the following line
//            String sample = "Sample " + line.substring(0, index);
            String sample = line.substring(0, index);
            String subType = sampleToType.get(sample);
            fu.printLine(sample + line.substring(index) + "\t" + subType);
        }
        fu.close();
    }
    
    @Test
    public void classifyBasedOnSSPs() throws Exception {
//        String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatio_UGID.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormZScore_UGID.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormGlobalZScore_070111_UGID.txt";
        // GSE4922
//        String srcFileName = R3Constants.BREAST_DIR + "GSE4922FilteredOnSamplesZScore_091409_UGID.txt";
        // GSE1456
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_z_020911_UGID.txt";
        // GSE3143
//        String srcFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_MappedGenes_z_012111_UGID.txt";
        // GSE1992
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1992/GSE1992_Gene_Exp_z_020711_UGID.txt";
        // GSE18229
        String srcFileName = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Gene_Exp_z_080311_UGID.txt";
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(srcFileName);
        // Load the centroids
        String centroidsFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/5centroids.txt";
        fu.setInput(centroidsFileName);
        List<Centroid> centroids = new ArrayList<Centroid>();
        String line = fu.readLine();
        String[] names = line.split("\t");
        for (int i = 2; i < names.length; i++) {
            Centroid centroid = new Centroid();
            centroid.name = names[i];
            centroids.add(centroid);
        }
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            for (int i = 2; i < tokens.length; i++) {
                Centroid centroid = centroids.get(i - 2);
                centroid.addGene(gene);
                centroid.addValue(new Double(tokens[i]));
            }
        }
        fu.close();
        // Get the list of all samples in NEJM
        fu.setInput(srcFileName);
        line = fu.readLine();
        fu.close();
        names = line.split("\t");
        System.out.println("Sample\tMax_Corr\tSub_Type");
        for (int i = 1; i < names.length; i++) {
            String sample = names[i];
            List<Double> sampleValues = new ArrayList<Double>();
            List<Double> centroidValues = new ArrayList<Double>();
            double maxValue = Double.NEGATIVE_INFINITY;
            Centroid found = null;
            for (Centroid centroid : centroids) {
                for (String gene : centroid.geneNames) {
                    Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
                    if (sampleToValue == null)
                        continue;
                    Double sampleValue = sampleToValue.get(sample);
                    if (sampleValue == null)
                        continue;
                    sampleValues.add(sampleValue);
                    centroidValues.add(centroid.getValue(gene));
                }
                // To avoid any error
                if (sampleValues.size() == 0 || centroidValues.size() == 0)
                    throw new IllegalStateException("Empty in expression values!");
                Double corr = MathUtilities.calculatePearsonCorrelation(sampleValues,
                                                                        centroidValues);
                centroidValues.clear();
                sampleValues.clear();
                if (corr > maxValue) {
                    maxValue = corr;
                    found = centroid;
                }
            }
            if (maxValue >= 0.10)
                System.out.println(sample + "\t" + maxValue + "\t" + found.name);
            else
                System.out.println(sample + "\t" + maxValue + "\tNA");
        }
    }
    
    @Test
    public void generateGSEAClassFile() throws IOException {
        // Use this file to get the ordered sample names
        String expFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormGlobalZScore_070111.txt";
        fu.setInput(expFileName);
        String line = fu.readLine();
        fu.close();
        // Just need the first line
        List<String> samples = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++)
            samples.add(tokens[i]);
//        // Quick check
//        for (String sample : samples)
//            System.out.println(sample);
        System.out.println("Total samples: " + samples.size());
        // The sub-type source
        String clinFileName = R3Constants.BREAST_DIR + "Nejm_Full_Clinical_Table_With_SubTypes_070411.txt";
        fu.setInput(clinFileName);
        line = fu.readLine();
        Map<String, String> sampleToType = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            sampleToType.put(tokens[0], tokens[tokens.length - 1]);
        }
        fu.close();
        for (String sample : samples) {
            String type = sampleToType.get(sample);
            System.out.println(sample + "\t" + type);
        }
        // Generate class file: two types are exported. One for LumA and another for non-LumA.
        String targetFileName = R3Constants.BREAST_DIR + "Nejm_LumA_Type.cls";
        fu.setOutput(targetFileName);
        fu.printLine(samples.size() + " 2 1");
        fu.printLine("# LumA Non_LumA");
        StringBuilder builder = new StringBuilder();
        for (String sample : samples) {
            String type = sampleToType.get(sample);
            if (type.equals("LumA"))
                builder.append("0");
            else
                builder.append("1");
            builder.append(" ");
        }
        fu.printLine(builder.toString());
        fu.close();             
    }
    
    @Test
    public void generateUGIDMapForNejm() throws IOException {
        String mapFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/NejmGenesToCLID.txt";
//        String mapFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE4922GenesToCLID.txt";
//        String mapFileName = R3Constants.BREAST_DIR + "MolecularPortraits_SSP/GSE1456GenesToCLID.txt";
        Map<String, String> nameToUGID = generateNameToUGID(mapFileName);
//        Map<String, String> nameToUGID = loadGeneNameToUGID();
        // This map is used to make sure only one UGIG should be output out.
        Map<String, Integer> ugidToNumber = new HashMap<String, Integer>();
        for (String name : nameToUGID.keySet()) {
            String ugid = nameToUGID.get(name);
            Integer number = ugidToNumber.get(ugid);
            if (number == null)
                ugidToNumber.put(ugid, 1);
            else
                ugidToNumber.put(ugid, ++number);
        }
        System.out.println("Total size: " + nameToUGID.size());
        Set<String> values = new HashSet<String>(nameToUGID.values());
        System.out.println("Total values: " + values.size());
        //String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatio.txt";
        //String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormZScore.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormGlobalZScore_070111.txt";
        String srcFileName = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Gene_Exp_z_080311.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "GSE4922FilteredOnSamplesZScore_091409.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_z_020911.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_MappedGenes_z_012111.txt";
//        String srcFileName = R3Constants.BREAST_DIR + "GSE1992/GSE1992_Gene_Exp_z_020711.txt";
        fu.setInput(srcFileName);
//        String targetFileName = R3Constants.BREAST_DIR + "NejmLogRatio_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormZScore_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "NejmLogRatioNormGlobalZScore_070111_UGID.txt";
        String targetFileName = R3Constants.BREAST_DIR + "GSE18229/GSE18229_Gene_Exp_z_080311_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE4922FilteredOnSamplesZScore_091409_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE1456/GSE1456_z_020911_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE3143/GSE3143_MappedGenes_z_012111_UGID.txt";
//        String targetFileName = R3Constants.BREAST_DIR + "GSE1992/GSE1992_Gene_Exp_z_020711_UGID.txt";
        fu.setOutput(targetFileName);
        String line = fu.readLine();
        fu.printLine(line);
        int length = line.split("\t").length;
        while ((line = fu.readLine()) != null) {
            int index = line.indexOf("\t");
            String geneName = line.substring(0, index);
            String ugid = nameToUGID.get(geneName);
            if (ugid == null)
                continue;
            Integer number = ugidToNumber.get(ugid);
            if (number > 1)
                continue; // Don't want to output duplicated UGIDs
            fu.printLine(ugid + line.substring(index));
        }
        fu.close();
    }

    private Map<String, String> generateNameToUGID(String mapFileName)
            throws IOException {
        Map<String, String> nameToUGID = new HashMap<String, String>();
        fu.setInput(mapFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 2)
                continue;
            if (tokens[1].length() == 0 || !tokens[1].startsWith("Hs."))
                continue;
            nameToUGID.put(tokens[0], tokens[1]);
        }
        fu.close();
        return nameToUGID;
    }
    
    private class Centroid {
        String name;
        List<String> geneNames;
        List<Double> values;
        
        public Centroid() {
            
        }
        
        public void addGene(String gene) {
            if (geneNames == null)
                geneNames = new ArrayList<String>();
            geneNames.add(gene);
        }
        
        public Double getValue(String gene) {
            int index = geneNames.indexOf(gene);
            return values.get(index);
        }
        
        public void addValue(Double value) {
            if (values == null)
                values = new ArrayList<Double>();
            values.add(value);
        }
    }
    
}
