/*
 * Created on Jan 25, 2011
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

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.r3.EntrezGeneAnalyzer;
import org.reactome.r3.util.GeneExpressionDataSet;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to handle downloaded soft text files from the GEO web site.
 * @author wgm
 *
 */
public class GSESoftTextDataHandler extends GSEDataHandler {
    // For control information
    // Control name in annotation
    private String controlColName = "CONTROL_TYPE"; // Default
    // Actual control name
    private String controlTypeName = "control";
    
    public GSESoftTextDataHandler() {
    }
    
    
    
    public String getControlColName() {
        return controlColName;
    }



    public void setControlColName(String controlColName) {
        this.controlColName = controlColName;
    }



    public String getControlTypeName() {
        return controlTypeName;
    }



    public void setControlTypeName(String controlTypeName) {
        this.controlTypeName = controlTypeName;
    }

    public Map<String, String> loadProbeIdToGene(String fileName, List<String> escapedPlatforms) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, String> probesetToGene = new HashMap<String, String>();
        boolean shouldRead = false;
        // Need to get the index for gene symbol
        int geneSymbolIndex = -1;
        int controlTypeIndex = -1;
        String platformId = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^PLATFORM")) {
                int index = line.indexOf("=");
                platformId = line.substring(index + 1).trim();
                geneSymbolIndex = -1;
                controlTypeIndex = -1;
            }
            else if (line.startsWith("!platform_table_begin")) {
                if (escapedPlatforms != null && escapedPlatforms.contains(platformId))
                    continue;
                // Should check next line
                line = fu.readLine();
//                System.out.println(line);
                String[] tokens = line.split("\t");
                for (int i = 0; i < tokens.length; i++) {
                    String token = tokens[i];
                    if (token.toLowerCase().equals("gene_symbol") ||
                            token.toLowerCase().equals("gene symbol")) {
                        geneSymbolIndex = i;
                    }
                    else if (token.equals(controlColName)) // Used to filter out controls
                        controlTypeIndex = i;
                }
                if (geneSymbolIndex == -1) {
                    throw new IllegalArgumentException("Cannot find GENE_SYMBOL column in the annotation section: " + platformId);
                }
                if (controlTypeIndex == -1) {
                    throw new IllegalArgumentException("Cannot find CONTROL_TYPE column in the annotation section: " + platformId);
                }
                //                System.out.println("Gene symbol index: " + geneSymbolIndex);
                shouldRead = true;
            }
            else if (line.startsWith("!platform_table_end")) {
                shouldRead = false;
                // Don't break the loop in case there are multiple platforms in the soft text file
                //break;
                // In case GPL1390, some special process is needed since no gene symbol has been provided.
                // In the original file, GB_ACC has been pretended as "GENE_SYMBOL" so that values in this
                // column can be loaded
                if (platformId.equals("GPL1390")) {
                    mapGPL1390GeneSymbols(probesetToGene);
                }
            }
            else if (shouldRead) {
                String[] tokens = line.split("\t");
                //                System.out.println(line);
                // Check CONTROL_TYPE if it is a control
                if (tokens.length > controlTypeIndex &&
                    tokens[controlTypeIndex].length() > 0 && 
                    tokens[controlTypeIndex].equals(controlTypeName)) // ||      
//                     !tokens[controlTypeIndex].equals("FALSE")))
                    continue;
                // Make sure there is GENE_SYMBOL in the line
                if (tokens.length > geneSymbolIndex && 
                        tokens[geneSymbolIndex].length() > 0 && // Not empty names
                        !tokens[geneSymbolIndex].contains("///") &&
                        !tokens[geneSymbolIndex].contains("|")) { // Don't want to probesets that can be mapped to multiple genes
                    probesetToGene.put(platformId + "." + tokens[0], tokens[geneSymbolIndex]);
                }
            }
        }
        fu.close();
        return probesetToGene;
    }

    /**
     * This method is used to load the mapping from probesets to gene symbols from a GEO soft text
     * file.
     * @param fileName
     * @return keys in the the returned map should be platform.probesetId, and values are gene symbols.
     * @throws IOException
     */
    @Override
    public Map<String, String> loadProbeIdToGene(String fileName) throws IOException {
        return loadProbeIdToGene(fileName, null);
    }
    
    private void mapGPL1390GeneSymbols(Map<String, String> platformProbeIdToSymbol) throws IOException {
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        Map<String, String> rnaAccToGeneSymbol = entrezAnalyzer.loadRNAAccToGeneSymbol();
        for (Iterator<String> it = platformProbeIdToSymbol.keySet().iterator(); it.hasNext();) {
            String platformProbeId = it.next();
            if (platformProbeId.startsWith("GPL1390")) {
                String acc = platformProbeIdToSymbol.get(platformProbeId);
                String symbol = rnaAccToGeneSymbol.get(acc);
                if (symbol == null)
                    it.remove();
                else
                    platformProbeIdToSymbol.put(platformProbeId, symbol);
            }
        }
    }
    
    /**
     * This method is used to extract a subset from an original loaded platformProbeIdToGene map.
     * @param platformProbeIdToGene
     * @return
     */
    private Map<String, String> getProbeIdToGene(Map<String, String> platformProbeIdToGene,
                                                 String platform) {
        Map<String, String> probeIdToGene = new HashMap<String, String>();
        int index = 0;
        for (String platformProbeId : platformProbeIdToGene.keySet()) {
            if (platformProbeId.startsWith(platform)) {
                index = platformProbeId.indexOf(".");
                probeIdToGene.put(platformProbeId.substring(index + 1),
                                  platformProbeIdToGene.get(platformProbeId));
            }
        }
//        for (String probeId : probeIdToGene.keySet()) {
//            String gene = probeIdToGene.get(probeId);
//            System.out.println(probeId + "\t" + gene);
//        }
        return probeIdToGene;
    }
    
    public Map<String, String> loadSampleFeature(String fileName,
                                                 String featureName,
                                                 String delimit) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, String> sampleToFeature = new HashMap<String, String>();
        String sample = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^SAMPLE")) { // Sample header
                int index = line.indexOf("=");
                sample = line.substring(index + 1).trim();
            }
            else if (line.startsWith(featureName)) {
                int index = line.indexOf(delimit);
                String feature = line.substring(index + 1).trim();
                if (!feature.equals("--")) { // Used for empty field
                    sampleToFeature.put(sample, feature);
                }
            }
        }
        fu.close();
        return sampleToFeature;
    }
    
    /**
     * This method is used to load all sample features.
     * @param fileName
     * @param featureTitle
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, String>> loadSampleToFeatures(String fileName,
                                                                 String featureTitle,
                                                                 String delimit) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, Map<String, String>> sampleToFeatureNameToValue = new HashMap<String, Map<String,String>>();
        String sample = null;
        Map<String, String> featureToValue = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^SAMPLE")) { // Sample header
                int index = line.indexOf("=");
                sample = line.substring(index + 1).trim();
                featureToValue = new HashMap<String, String>();
                sampleToFeatureNameToValue.put(sample, featureToValue);
            }
            else if (line.startsWith(featureTitle)) {
                int index = line.indexOf(delimit);
                String featureValue = line.substring(index + 1).trim();
                index = featureValue.indexOf(":");
                String feature = featureValue.substring(0, index).trim();
                String value = featureValue.substring(index + 1).trim();
                if (value.equals("--"))
                    value = null;
                featureToValue.put(feature, value); // Save a null value too.
            }
        }
        fu.close();
        return sampleToFeatureNameToValue;
    }
    
    /**
     * Note: This is very important!!! The following loading code assuming all features are sorted in the same order. 
     * Otherwise, the results will not be right.
     * @param fileName
     * @param platformId
     * @return
     * @throws IOException
     */
    public GeneExpressionDataSet loadGeneExpressionData(String fileName,
                                                        String platformId) throws IOException {
        GeneExpressionDataSet dataset = new GeneExpressionDataSet();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        boolean isInData = false;
        List<List<Double>> allValues = null;
        boolean shouldRead = false;
        String sample = null;
        // Use a map to load feature to value pair in case the order of features have been changed.
        // This is very important!!!
        Map<String, Double> featureToValue = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^SAMPLE")) {
                index = line.indexOf("=");
                sample = line.substring(index + 1).trim();
            }
            else if (line.startsWith("!Sample_platform_id")) {
                index = line.indexOf("=");
                String tmpPlatformId = line.substring(index + 1).trim();
                if (tmpPlatformId.equals(platformId)) {
                    shouldRead = true;
                    // Default to here to add samples
                    dataset.addSample(sample);
                    // Reset these flags
                    isInData = false;
                }
            }
            else if (shouldRead && line.equals("!sample_table_begin")) {
                // Escape the title line
                line = fu.readLine();
                isInData = true;
                featureToValue.clear();
            }
            else if (shouldRead && line.equals("!sample_table_end")) {
                if (dataset.getFeatureList() == null) {
                    List<String> featureList = new ArrayList<String>(featureToValue.keySet());
                    Collections.sort(featureList);
                    dataset.setFeatureList(featureList);
                }
                else {
                    // make sure the feature list is the same
                    List<String> currentFeatures = new ArrayList<String>(featureToValue.keySet());
                    Collections.sort(currentFeatures);
                    if (!dataset.getFeatureList().equals(currentFeatures)) {
                        throw new IllegalStateException(sample + " has a different feature order!");
                    }
                    System.out.println(sample + " has the same feature order!");
                }
                // Dump values
                if (allValues == null) {
                    allValues = new ArrayList<List<Double>>();
                    for (int i = 0; i < dataset.getFeatureList().size(); i++)
                        allValues.add(new ArrayList<Double>());
                }
                List<String> features = dataset.getFeatureList();
                for (int i = 0; i < features.size(); i++) {
                    List<Double> featureValues = allValues.get(i);
                    String feature = features.get(i);
                    Double value = featureToValue.get(feature);
                    featureValues.add(value);
                }
                isInData = false;
                shouldRead = false;
//                for (int i = 0; i < 10; i++)
//                    System.out.println("######");
            }
            else if (isInData) {
//                System.out.println(line);
                String[] tokens = line.split("\t");
                if (tokens[0].startsWith("AFFX-"))
                    continue; // Affy controls. Just ignore them!
                if (tokens[1].length() == 0) 
                    featureToValue.put(tokens[0], null);
                else
                    featureToValue.put(tokens[0], new Double(tokens[1]));
            }
        }
        dataset.setValues(allValues);
        return dataset;
    }
    
    @Test
    public void testLoadExpressionDataset() throws IOException {
        String fileName = "datasets/BreastCancer/GSE1922/GSE1992_family.soft";
        GeneExpressionDataSet dataset = loadGeneExpressionData(fileName, "GPL887");
        System.out.println("Total samples: " + dataset.getSampleList().size());
        System.out.println("Total features: " + dataset.getFeatureList().size());
        // Just a quick check
        List<List<Double>> values = dataset.getValues();
        List<Double> featureValues = values.get(0);
        for (Double v : featureValues)
            System.out.println(v);
    }

    @Test
    public void testLoadProbesetToGenes() throws IOException {
        //String fileName = "datasets/BreastCancer/GSE1922/GSE1992_family.soft";
        String fileName = "datasets/BreastCancer/GSE18229/GSE18229_family.soft";
        Map<String, String> idToGene = loadProbeIdToGene(fileName);
        System.out.println("Total probeset: " + idToGene.size());
        int count = 0;
        for (String id : idToGene.keySet()) {
            String gene = idToGene.get(id);
            System.out.println(id + "\t" + gene);
            count ++;
            if (count == 100)
                break;
        }
    }
    
    @Test
    public void checkSoftFile() throws IOException {
//        String fileName = "datasets/BreastCancer/GSE1992/GSE1992_family.soft";
//        fileName = CancerConstants.BREAST_DIR + "GSE1456/GSE1456_family.soft";
        String fileName = CancerConstants.BREAST_DIR + "GSE18229/GSE18229_family.soft";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^"))
                System.out.println(line);
//            if (line.startsWith("^PLATFORM"))
//                System.out.println(line);
        }
        fu.close();
        Map<String, String> probeIdToGene = loadProbeIdToGene(fileName);
        System.out.println("Total probeids: " + probeIdToGene.size());
        // Recapture the list of platformts
        Set<String> platforms = new HashSet<String>();
        int index = 0;
        for (String probeId : probeIdToGene.keySet()) {
            index = probeId.indexOf(".");
            platforms.add(probeId.substring(0, index));
        }
        System.out.println("Platforms: " + platforms);
//        GeneExpressionDataSet dataset = loadGeneExpressionData(fileName, "GPL1390");
//        System.out.println("GPL1390 probeids: " + dataset.getFeatureList().size());
    }
    
    /**
     * This helper method is used to get the list of platforms used in a soft file.
     * @param fileName
     * @return
     * @throws IOException
     */
    public List<String> loadPlatforms(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = null;
        List<String> platforms = new ArrayList<String>();
        int index = 0;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("^PLATFORM")) {
                index = line.indexOf("=");
                String platform = line.substring(index + 1).trim();
                platforms.add(platform);
            }
        }
        fu.close();
        return platforms;
    }
    
    @Test
    public void testLoadSampleToFeature() throws IOException {
        String dirName = "datasets/BreastCancer/GSE1922/";
        String fileName = "datasets/BreastCancer/GSE1922/GSE1992_family.soft";
        String featureName = "!Sample_characteristics_ch2 = overall survival event (0=alive, 1=dod or doc)";
        Map<String, String> sampleToEvent = loadSampleFeature(fileName, featureName, ":");
        featureName = "!Sample_characteristics_ch2 = overall suvival months";
        Map<String, String> sampleToTime = loadSampleFeature(fileName, featureName, ":");
        List<String> sampleList = new ArrayList<String>(sampleToTime.keySet());
        Collections.sort(sampleList);
        String output = dirName + "GSE1922_OS_Info.txt";
        fu.setOutput(output);
        fu.printLine("Sample_ID\tOSEVENT\tOSDURATION");
        for (String sample : sampleList) {
            String event = sampleToEvent.get(sample);
            String time = sampleToTime.get(sample);
            System.out.println(sample + "\t" + event + "\t" + time);
            fu.printLine(sample + "\t" + event + "\t" + time);
        }
        fu.close();
    }
    
    @Test
    public void processGSE1992SoftTextFile() throws IOException {
        String dirName = "datasets/BreastCancer/GSE1992/";
        String fileName = dirName + "GSE1992_family_changed.soft";
        //String outFileName = dirName + "GSE1992_Gene_Exp_020711.txt";
//        String outFileName = dirName + "GSE1992_Gene_Exp_z_080311.txt";
//        String outFileName = dirName + "GSE1992_Gene_Exp_080311.txt";
//        String outFileName = dirName + "GSE1992_GPL1390_Gene_Exp_080411.txt";
//        String outFileName = dirName + "GSE1992_GPL885_887_Gene_Exp_080411.txt";
        String outFileName = dirName + "GSE1992_GPL1390_z_Gene_Exp_080411.txt";
        Map<String, String> platformProbeIdToGene = loadProbeIdToGene(fileName);
        List<String> platforms = loadPlatforms(fileName);
        //platforms.remove("GPL1390");
        platforms.clear();
        platforms.add("GPL1390");
        // For this data set, samples in different platforms are different.
        // Have to merge samples in different platforms manually
        processSoftTextFile(fileName, 
                            outFileName, 
                            platformProbeIdToGene,
                            platforms,
                            false);
    }
    
    @Test
    public void processGSE9899SoftTextFile() throws IOException {
        controlColName = "Sequence Type";
        controlTypeName = "Control sequence";
        String dirName = "datasets/TCGA/OvarianCancer/TothillDataset/";
//        String fileName = dirName + "GSE9899_family.soft";
        String fileName = dirName + "GSE9891_family.soft";
//        String outFileName = dirName + "GSE9899_z_Gene_Exp_111811.txt";
//        String outFileName = dirName + "GSE9899_Gene_Exp_111811_no_avg.txt";
//        String outFileName = dirName + "GSE9891_z_Gene_Exp_111811.txt";
//        String outFileName = dirName + "GSE9891_z_Gene_Exp_111811.txt";
        String outFileName = dirName + "GSE9891_Gene_Exp_111811.txt";
        Map<String, String> platformProbeIdToGene = loadProbeIdToGene(fileName);
        String platform = "GPL570";
        String[] testIds = new String[] {
                platform + ".1569555_at",
                platform + ".224209_s_at"
        };
        for (String testId : testIds)
            System.out.println(testId + "\t" + platformProbeIdToGene.get(testId));
        List<String> platforms = loadPlatforms(fileName);
        System.out.println("Total platforms: " + platforms.size());
        processSoftTextFile(fileName, 
                            outFileName, 
                            platformProbeIdToGene,
                            platforms,
                            false);
    }
    
    @Test
    public void processGSE26712SoftTextFile() throws IOException {
        controlColName = "Sequence Type";
        controlTypeName = "Control sequence";
        String dirName = "datasets/TCGA/OvarianCancer/GSE26712/";
        String fileName = dirName + "GSE26712_family.soft";
        String outFileName = dirName + "GSE26712_Gene_Exp_111911.txt";
        Map<String, String> platformProbeIdToGene = loadProbeIdToGene(fileName);
        List<String> platforms = loadPlatforms(fileName);
        System.out.println("Total platforms: " + platforms.size());
        processSoftTextFile(fileName, 
                            outFileName, 
                            platformProbeIdToGene,
                            platforms,
                            false);
    }
    
    /**
     * Do a zscore transformation without considering 10 normal samples in the GSE26712 dataset.
     * @throws IOException
     */
    @Test
    public void ztransformGSE26712ExpresionNoNormalSamples() throws IOException {
        // GSE657519 - GSE657528 are normal tissue controls
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        String dirName = "datasets/TCGA/OvarianCancer/GSE26712/";
        String srcFileName = dirName + "GSE26712_Gene_Exp_111911.txt";
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(srcFileName);
        System.out.println("Total genes: " + geneToSampleToValue.size());
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            for (String sample : sampleToValue.keySet()) {
                // Do a quick check about name
                int number = Integer.parseInt(sample.substring(3));
                if (number >= 657519 && number <= 657528)
                    continue;
                stat.addValue(sampleToValue.get(sample));
            }
        }
        System.out.println("SD: " + stat.getStandardDeviation());
        System.out.println("Mean: " + stat.getMean());
//        if (true)
//            return;
        double sd = stat.getStandardDeviation();
        double mean = stat.getMean();
        // Reset the value
        for (String gene : geneToSampleToValue.keySet()) {
            Map<String, Double> sampleToValue = geneToSampleToValue.get(gene);
            for (String sample : sampleToValue.keySet()) {
                // Do a quick check about name
                int number = Integer.parseInt(sample.substring(3));
                if (number >= 657519 && number <= 657528)
                    continue;
                Double value = sampleToValue.get(sample);
                value = (value - mean) / sd;
                sampleToValue.put(sample, value);
            }
        }
//        if (true)
//            return;
        helper.outputGeneExp(geneToSampleToValue, dirName + "GSE26712_no_normal_z_Gene_Exp_112111.txt");
    }
    
    @Test
    public void processGSE26712ClinFile() throws IOException {
        String dirName = "datasets/TCGA/OvarianCancer/GSE26712/";
        String fileName = dirName + "GSE26712_family.soft";
        Map<String, String> sampleToStatus = loadSampleFeature(fileName, 
                                                              "!Sample_characteristics_ch1 = status", 
                                                              ": ");
        Map<String, String> sampleToSurvival = loadSampleFeature(fileName,
                                                                 "!Sample_characteristics_ch1 = survival years",
                                                                 ": ");
        System.out.println("Sample\tOSDURATION_YEARS\tStatus\tOSEVENT");
        for (String name : sampleToStatus.keySet()) {
            String status = sampleToStatus.get(name);
            String survival = sampleToSurvival.get(name);
            if (status == null && survival == null)
                continue;
            if (status == null)
                status = "";
            if (survival == null)
                survival = "";
            String death = null;
            if (status.startsWith("DOD"))
                death = "1";
            else
                death = "0";
            System.out.println(name + "\t" + survival + "\t" + status + "\t" + death);
        }
    }
    
    private void checkDataset(GeneExpressionDataSet dataset) {
        List<String> samples = dataset.getSampleList();
        System.out.println("First sample: " + samples.get(0));
        List<String> features = dataset.getFeatureList();
        for (int i = 0; i < features.size(); i++) {
            List<Double> values = dataset.getValues().get(i);
            System.out.println(features.get(i) + "\t" + values.get(0));
        }
    }

    protected void processSoftTextFile(String fileName,
                                       String outFileName,
                                       Map<String, String> platformProbeIdToGene,
                                       List<String> platforms,
                                       boolean needZscoreTransform) throws IOException {
        List<GeneExpressionDataSet> datasets = new ArrayList<GeneExpressionDataSet>();
        for (String platform : platforms) {
            System.out.println("Platform: " + platform);
            GeneExpressionDataSet dataset = loadGeneExpressionData(fileName, platform);
//            checkDataset(dataset);
            System.out.println("Total samples: " + dataset.getSampleList().size());
            Map<String, String> probeIdToGene = getProbeIdToGene(platformProbeIdToGene,
                                                                 platform);
//            String[] testIds = new String[] {
//                    "1569555_at",
//                    "224209_s_at"
//            };
//            for (String id : testIds)
//                System.out.println(id + "\t" + probeIdToGene.get(id));
            mapProbesetToGenes(dataset, probeIdToGene);
            dataset = averageValuesForSameGenes(dataset);
            datasets.add(dataset);
        }
        GeneExpressionDataSet mergedDataset = mergeDatasetsWithDiffSamples(datasets);
        if (needZscoreTransform)
            mergedDataset.zscoreTansformation(); // Do a global z-score transformation
//        mergedDataset.geneWiseMedianTransformation();
//        mergedDataset.sampleWiseZscoreTransformation();
        mergedDataset.export(outFileName);
        System.out.println("Total samples: " + mergedDataset.getSampleList().size());
        System.out.println("Total genes: " + mergedDataset.getFeatureList().size());
    }
    
    @Test
    public void processGSE18229SoftTextFile() throws IOException {
        String dirName = CancerConstants.BREAST_DIR + "GSE18229/";
        String fileName = dirName + "GSE18229_family.soft";
        String platformFileName = dirName + "GSE18229_family_platforms.soft";
//        String outFileName = dirName + "GSE18229_Gene_Exp_z_080311.txt";
//        String outFileName = dirName + "GSE18229_GPL1390_Gene_Exp_080311.txt";
        String outFileName = dirName + "GSE18229_Gene_Exp_111811.txt";
        Map<String, String> platformProbeIdToGene = loadProbeIdToGene(platformFileName);
        List<String> platforms = loadPlatforms(platformFileName);
        System.out.println("Total platforms: " + platforms.size());
//        platforms.clear();
//        platforms.add("GPL1390");
        processSoftTextFile(fileName, 
                            outFileName, 
                            platformProbeIdToGene, 
                            platforms,
                            true);
    }
    
    public GeneExpressionDataSet mergeDatasetsWithDiffSamples(List<GeneExpressionDataSet> datasets) {
        // Convert GeneExpressionDataSet into Map of samples to gene to expression values
        List<Map<String, Map<String, Double>>> sampleToGeneToValueList = new ArrayList<Map<String, Map<String, Double>>>();
        for (GeneExpressionDataSet dataset : datasets) {
            Map<String, Map<String, Double>> sampleToGeneToValue = dataset.convertToSampleToFeatureToValueMap();
            sampleToGeneToValueList.add(sampleToGeneToValue);
        }
        // Merge three maps. Since all maps have different keys (samples here), we can just merge them
        // together without further process
        Map<String, Map<String, Double>> mergedMap = new HashMap<String, Map<String,Double>>();
        for (Map<String, Map<String, Double>> sampleToGeneToValue : sampleToGeneToValueList) {
            mergedMap.putAll(sampleToGeneToValue);
        }
        // Convert the Map into GeneExpressionDataSet to use its methods
        GeneExpressionDataSet mergedDataset = GeneExpressionDataSet.convertSampleToFeatureToValueMap(mergedMap);
        return mergedDataset;
    }
    
    @Test
    public void mapSamplesInClinInfoFile() throws IOException {
        String dirName = CancerConstants.BREAST_DIR + "GSE18229/";
        String fileName = dirName + "GSE18229_family.soft";
        Map<String, String> sampleToFeature = loadSampleFeature(fileName,
                                                                "!Sample_title",
                                                                "=");
        Map<String, String> titleToSample = new HashMap<String, String>();
        for (String sample : sampleToFeature.keySet()) {
            String feature = sampleToFeature.get(sample);
            titleToSample.put(feature, sample);
        }
        String clinInfoFileName = dirName + "GSE18229ClinicalInfo_Original.txt";
        String outClinFileName = dirName + "GSE18229_Clin_info.txt";
        fu.setInput(clinInfoFileName);
        fu.setOutput(outClinFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        for (int i = 2; i < tokens.length; i++)
            builder.append("\t").append(tokens[i]);
        fu.printLine(builder.toString());
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String gsmId = titleToSample.get(tokens[1]);
            if (gsmId == null) {
                System.out.println(tokens[1] + " has no map!");
            }
            builder.setLength(0);
            builder.append(gsmId);
            for (int i = 2; i < tokens.length; i++)
                builder.append("\t").append(tokens[i]);
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    @Test
    public void generateGSE1992FullSampleFeatureFile() throws IOException {
        String dirName = "datasets/BreastCancer/GSE1992/";
        String fileName = dirName + "GSE1992_family.soft";
        String outFileName = dirName + "GSE1992_Clin_Info.txt";
        Map<String, Map<String, String>> sampleToFeatureToValue = loadSampleToFeatures(fileName, 
                                                                                       "!Sample_characteristics_ch2", 
                                                                                       "=");
        // Output
        fu.setOutput(outFileName);
        // Generate a table header
        StringBuilder builder = new StringBuilder();
        builder.append("Sample");
        Set<String> features = new HashSet<String>();
        // Load the whole data set once in case any feature list inconsistence
        for (Map<String, String> featureToValue : sampleToFeatureToValue.values()) {
            features.addAll(featureToValue.keySet());
        }
        List<String> featureList = new ArrayList<String>(features);
        Collections.sort(featureList);
        for (String feature : featureList)
            builder.append("\t").append(feature);
        fu.printLine(builder.toString());
        builder.setLength(0);
        List<String> sampleList = new ArrayList<String>(sampleToFeatureToValue.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            Map<String, String> featureToValue = sampleToFeatureToValue.get(sample);
            builder.append(sample);
            for (String feature : featureList) {
                String value = featureToValue.get(feature);
                builder.append("\t");
                if (value == null)
                    builder.append("");
                else
                    builder.append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This method is used to process the download GSE1456_family.soft text files to generate
     * consolidated gene, sample, expression values. 
     * Note: there are two platforms bundled together in this single file.
     * @throws IOException
     */
    @Test
    public void processGSE1456SoftTextFile() throws IOException {
        String dirName = "datasets/BreastCancer/GSE1456/";
        String fileName = dirName + "GSE1456_family.soft";
        String outFileName = dirName + "GSE1456_z_020911.txt";
        GeneExpressionDataSet gpl96Dataset = loadGeneExpressionData(fileName, "GPL96");
        System.out.println("Total samples in GPL96: " + gpl96Dataset.getSampleList().size());
        System.out.println("Total features in GPL96: " + gpl96Dataset.getFeatureList().size());
        GeneExpressionDataSet gpl97Dataset = loadGeneExpressionData(fileName, "GPL97");
        System.out.println("Total samples in GPL97: " + gpl97Dataset.getSampleList().size());
        System.out.println("Total features in GPL97: " + gpl97Dataset.getFeatureList().size());
        Map<String, String> probeIdToGene = loadProbeIdToGene(fileName);
        System.out.println("Total probe ids: " + probeIdToGene.size());
        Map<String, String> gsmIdToSampleId = loadSampleFeature(fileName, 
                                                                "!Sample_title", 
                                                                "=");
        // Need some clean up to remove the postfix
        for (String gsmId: gsmIdToSampleId.keySet()) {
            String sampleId = gsmIdToSampleId.get(gsmId);
            int index = sampleId.indexOf(" ");
            sampleId = sampleId.substring(0, index).trim();
            gsmIdToSampleId.put(gsmId, sampleId);
        }
        // Reset the samples list in the two GeneExpressionDataSet objects
        mapSamples(gpl96Dataset, 
                   gsmIdToSampleId);
        mapSamples(gpl97Dataset,
                   gsmIdToSampleId);
        // Merge two GeneExpressionDataSet objects
        gpl96Dataset.merge(gpl97Dataset);
        mapProbesetToGenes(gpl96Dataset, probeIdToGene);
        GeneExpressionDataSet mergedDataSet = averageValuesForSameGenes(gpl96Dataset);
        mergedDataSet.zscoreTansformation(); // Do a global z-score transformation
        mergedDataSet.export(outFileName);
    }

    private void mapSamples(GeneExpressionDataSet gpl96Dataset,
                            Map<String, String> gsmIdToSampleId) {
        List<String> gpl96Samples = gpl96Dataset.getSampleList();
        List<String> mappedGPL96Samples = new ArrayList<String>(gpl96Samples.size());
        for (String sample : gpl96Samples) {
            String sampleId = gsmIdToSampleId.get(sample);
            mappedGPL96Samples.add(sampleId);
        }
        gpl96Dataset.setSampleList(mappedGPL96Samples);
    }
    
    /**
     * There are two platforms are bundled together in the GSE1456 soft file. This
     * method is used to do some sanity checkings.
     * @throws IOException
     */
    @Test
    public void checkGSE1456SoftFamilyFile() throws IOException {
        String fileName = "datasets/BreastCancer/GSE1456/GSE1456_family.soft";
        fu.setInput(fileName);
        String line = null;
        Set<String> probesets1 = new HashSet<String>();
        Set<String> probesets2 = new HashSet<String>();
        int index = 0;
        boolean shouldRead = false;
        Set<String> tmp = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("!platform_table_begin")) {
                shouldRead = true;
                if (tmp == null)
                    tmp = probesets1;
                else
                    tmp = probesets2;
            }
            else if (line.startsWith("!platform_table_end")) {
                shouldRead = false;
                // Don't break the loop in case there are multiple platforms in the soft text file
                //break;
            }
            else if (shouldRead) {
                String[] tokens = line.split("\t");
                tmp.add(tokens[0]);
            }
        }
        fu.close();
        System.out.println("Probesets in platform 1: " + probesets1.size());
        System.out.println("Probesets in platform 2: " + probesets2.size());
        Set<String> total = new HashSet<String>();
        total.addAll(probesets1);
        total.addAll(probesets2);
        System.out.println("Total: " + total.size());
        Set<String> shared = InteractionUtilities.getShared(probesets1, probesets2);
        System.out.println("Shared: " + shared.size());
        for (String id : shared)
            System.out.println(id);
        // Check sample titles
        Map<String, String> sampleToTitle = loadSampleFeature(fileName, "!Sample_title", "=");
        System.out.println("Total samples: " + sampleToTitle.size());
        List<String> sampleList = new ArrayList<String>(sampleToTitle.keySet());
        Collections.sort(sampleList);
        for (String sample : sampleList) {
            String title = sampleToTitle.get(sample);
            System.out.println(sample + "\t" + title);
        }
    }
}
