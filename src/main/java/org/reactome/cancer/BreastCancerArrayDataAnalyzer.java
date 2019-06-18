/*
 * Created on May 30, 2008
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.data.ProteinSequenceHandler;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.r3.EntrezGeneAnalyzer;
import org.reactome.r3.graph.GraphComponent;
import org.reactome.r3.graph.MutualInformationScoreCalculator;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This method is used to calculate z-scores for array data set.
 * @author wgm
 *
 */
public class BreastCancerArrayDataAnalyzer {
    private final String dirName = "datasets/BreastCancer/";
    private FileUtility fu = new FileUtility();
    
    public BreastCancerArrayDataAnalyzer() {
    }
    
    @Test
    public void clusterSamplesBasedOnFIs() throws Exception {
        String fileName = dirName + "GSE2034FilteredNormZScore.txt";
        // Calculate mutation information
        MutualInformationScoreCalculator miCalculator = new MutualInformationScoreCalculator();
        // Search components based on Wang et al data set
        miCalculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        miCalculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        List<String> samples = miCalculator.getSamples();
        Map<String, Set<String>> sampleToTopGenes = new HashMap<String, Set<String>>();
        final Map<String, List<Double>> geneToValues = miCalculator.getIdToValues();
        List<String> geneList = new ArrayList<String>(geneToValues.keySet());
        int slice = (int) (geneList.size() * 0.001);
        for (int i = 0; i < samples.size(); i++) {
            final int index = i;
            Collections.sort(geneList, new Comparator<String>() {
                public int compare(String gene1, String gene2) {
                    List<Double> list1 = geneToValues.get(gene1);
                    List<Double> list2 = geneToValues.get(gene2);
                    Double value1 = list1.get(index);
                    Double value2 = list2.get(index);
                    return value2.compareTo(value1);
                }
            });
            Set<String> topGenes = new HashSet<String>();
            // Pick the first 0.5% and last 0.5%
            for (int j = 0; j < slice; j++)
                topGenes.add(geneList.get(j));
            for (int j = geneList.size() - 1; j < geneList.size() - slice - 1; j --)
                topGenes.add(geneList.get(j));
            sampleToTopGenes.put(samples.get(i), topGenes);
        }
        // Try to do a network clustering
        NatureGBMAnalyzer gbmAnalyzer = new NatureGBMAnalyzer();
        gbmAnalyzer.hierchicallyClusterSamples(new ArrayList<String>(sampleToTopGenes.keySet()), 
                                               sampleToTopGenes);
    }
    
    /**
     * This method is used to generate an arff file for weka to do classifying.
     * @throws Exception
     */
    @Test
    public void generateArffFileForGenes() throws Exception {
        int geneNumber = 30;
        MutualInformationScoreCalculator calculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        calculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        calculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        String fileName = dirName + "GSE2034Genes30.arff";
//        // GSE4922
//        calculator.setSampleToPhenotypeInfo(loadGSE4922PatientInfo());
//        calculator.initIdToValueInfo(dirName + "GSE4922FilteredOnSamplesZScore.txt");
//        String fileName = dirName + "Components2_1_GSE4922_NO_P_I_cutoff_10_GSE4922.arff";        
//        // For van de Vijver data set
//        calculator.setSampleToPhenotypeInfo(dataAnalyzer.loadNejmPatientInfo());
//        calculator.initIdToValueInfo(dirName + "NejmLogRatioNormZScore.txt");
//        String fileName = dirName + "Components2_1_Nejm_NO_P_I_cutoff_11_Nejm.arff";

        final Map<String, Double> id2Mi = getMIForGene();
        List<String> sortedGenes = new ArrayList<String>(id2Mi.keySet());
        Collections.sort(sortedGenes, new Comparator<String>() {
            public int compare(String id1, String id2) {
                Double mi1 = id2Mi.get(id1);
                Double mi2 = id2Mi.get(id2);
                return mi2.compareTo(mi1);
            }
        });
        for (String gene : sortedGenes) {
            System.out.println(gene + "\t" + id2Mi.get(gene));
        }
        if (true)
            return;
        List<String> samples = calculator.getSamples();
        Map<String, String> sampleToType = calculator.getSampleToType();
        Map<String, List<Double>> idToValues = calculator.getIdToValues();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("@relation BreastCancerType\n\n");
        builder.append("@attribute CellType {true, false}\n");
        for (int i = 0; i < geneNumber; i++) {
            builder.append("@attribute ").append(sortedGenes.get(i));
            builder.append(" NUMERIC\n");
        }
        builder.append("\n@data");
        fu.printLine(builder.toString());
        builder.setLength(0);
        // Output the table
        for (int i = 0; i < samples.size(); i++) {
            String sample = samples.get(i);
            String type = sampleToType.get(sample);
            if (type.equals("0"))
                type = "false";
            else
                type = "true";
            builder.append(type);
            for (int j = 0; j < geneNumber; j++) {
                String gene = sortedGenes.get(j);
                List<Double> values = idToValues.get(gene);
                builder.append(",").append(values.get(i));
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * Used to check individual gene MI scores
     * @throws Exception
     */
    private Map<String, Double> getMIForGene() throws Exception {
        MutualInformationScoreCalculator calculator = new MutualInformationScoreCalculator();
        // For Wang et al data set
        calculator.initSampleToPhenotypeInfo(dirName + "GSE2034_Patient_Info.txt");
        calculator.initIdToValueInfo(dirName + "GSE2034FilteredNormZScore.txt");
        Map<String, List<Double>> idToValues = calculator.getIdToValues();
        Map<String, Double> id2mi = new HashMap<String, Double>();
        for (String id : idToValues.keySet()) {
            GraphComponent comp = new GraphComponent();
            comp.addNode(id);
            double score = calculator.calculateScore(comp);
            id2mi.put(id, score);
        }
        return id2mi;
    }
    
    /**
     * This method is used to calculate activities of clusters.
     * @throws IOException
     */
    @Test
    public void calculateActivitiesOfClusters() throws Exception {
        Map<Integer, Set<String>> clusterToGenes = getClusterToGenes();
        String fileName = dirName + "GSE2034FilteredZScore.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> samples = new ArrayList<String>();
        extractGSE2034Samples(line, samples);
        // Need to load gene expression values
        Map<String, List<String>> geneToValues = new HashMap<String, List<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            List<String> values = new ArrayList<String>();
            for (int i = 1; i < tokens.length; i++) 
                values.add(tokens[i]);
            geneToValues.put(tokens[0],
                             values);
        }
        fu.close();
        // Calculate activities
        Map<Integer, List<String>> clusterToValues = new HashMap<Integer, List<String>>();
        for (Iterator<Integer> it = clusterToGenes.keySet().iterator(); it.hasNext();) {
            Integer cluster = it.next();
            Set<String> genesInCluster = clusterToGenes.get(cluster);
            // Have to make sure more than half of genes having expression values
            int c = 0;
            for (String gene : genesInCluster) {
                if (geneToValues.containsKey(gene))
                    c ++;
            }
            if (((double)c / genesInCluster.size()) < 0.5d)
                continue; // Don't want to use it
            List<String> clusterAverage = averageValues(new ArrayList<String>(genesInCluster),
                                                        geneToValues);
            clusterToValues.put(cluster,
                                clusterAverage);
        }
        // Output
        fileName = dirName + "GSE2034ClusterValues.txt";
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Cluster");
        for (String sample : samples)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        for (Iterator<Integer> it = clusterToValues.keySet().iterator(); it.hasNext();) {
            Integer cluster = it.next();
            List<String> values = clusterToValues.get(cluster);
            builder.setLength(0);
            builder.append(cluster);
            for (String value : values)
                builder.append("\t").append(value);
            fu.printLine(builder.toString());
        }
        fu.close();
    }
        
    /**
     * This method is used to analyze the mapping from probeset to protein UniProt ids.
     * @throws Exception
     */
    @Test
    public void checkProbesetToUniprotMap() throws Exception {
        String fileName = dirName + "HG-U133A.na25.annot.csv";
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Set<String>> probesetToUniprotIds = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\",");
            String probeset = removeQuotationMarks(tokens[0]);
            String idString = removeQuotationMarks(tokens[19]);
            tokens = idString.split("///");
            Set<String> ids = new HashSet<String>();
            for (String id : tokens)
                ids.add(id.trim());
            probesetToUniprotIds.put(probeset, ids);
        }
        System.out.println("Total probeset: " + probesetToUniprotIds.size());
        checkSingleMappedUniProtIds(probesetToUniprotIds);
        // Try to consolidate UniProt ids
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        Map<String, Set<String>> filterProbeToProtIds = new HashMap<String, Set<String>>();
        for (Iterator<String> it = probesetToUniprotIds.keySet().iterator(); it.hasNext();) {
            String probeset = it.next();
            Set<String> ids = probesetToUniprotIds.get(probeset);
            if (ids.size() == 1) {
                filterProbeToProtIds.put(probeset, ids);
            }
            else {
                Set<String> filteredIds = seqHandler.consolidateProteinIds(ids);
                filterProbeToProtIds.put(probeset,
                                         filteredIds);
            }
        }
        checkSingleMappedUniProtIds(filterProbeToProtIds);
    }
    
    private void checkSingleMappedUniProtIds(Map<String, Set<String>> probesetToUniprotIds) {
        int count = 0;
        for (Iterator<String> it = probesetToUniprotIds.keySet().iterator(); it.hasNext();) {
            String probeset = it.next();
            Set<String> ids = probesetToUniprotIds.get(probeset);
            if (ids.size() == 1)
                count ++;
        }
        System.out.println("Single mapped: " + count);
    }
    
    /**
     * This method is used to check shared genes between two data sets.
     * @throws IOException
     */
    @Test
    public void checkSharedGenes() throws IOException {
        Set<String> nejmGeneSet = getNejmGenes();
        Set<String> gpl96GeneSet = getGPL96Genes();
        System.out.println("Nejm Genes: " + nejmGeneSet.size());
        System.out.println("GPL 96 genes: " + gpl96GeneSet.size());
        Set<String> copy = new HashSet<String>(gpl96GeneSet);
        copy.retainAll(nejmGeneSet);
        System.out.println("Shared genes: " + copy.size());
        // After the normalization
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        nejmGeneSet = entrezAnalyzer.normalizeGeneNames(nejmGeneSet);
        gpl96GeneSet= entrezAnalyzer.normalizeGeneNames(gpl96GeneSet);
        System.out.println("After entrez gene normalization:");
        System.out.println("Nejm genes: " + nejmGeneSet.size());
        System.out.println("GPL 96 genes: " + gpl96GeneSet.size());
        copy = new HashSet<String>(gpl96GeneSet);
        copy.retainAll(nejmGeneSet);
        System.out.println("Shared genes: " + copy.size());
        String intFileName = R3Constants.RESULT_DIR + "FI73InGene_061008.txt";
        Set<String> fis = fu.loadInteractions(intFileName);
        Set<String> geneInFis = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes in network: " + geneInFis.size());
        geneInFis.retainAll(gpl96GeneSet);
        System.out.println("Total genes shared between network and GSE2034: " + geneInFis.size());
        geneInFis.retainAll(nejmGeneSet);
        System.out.println("Total genes shared in network, GSE2034 and nejm: " + geneInFis.size());
    }
    
    private Set<String> getGPL96Genes() throws IOException {
        String fileName = dirName + "HG-U133A.na25.annot.csv";
        Set<String> geneNames = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\",");
            if (tokens.length < 15)
                continue;
            String nameString = removeQuotationMarks(tokens[14]);
            if (nameString.length() == 0 ||
                nameString.contains("///"))
                continue; // Multiple mapping
            geneNames.add(nameString.toUpperCase());
        }
        fu.close();
        return geneNames;
    }
    
    /**
     * This method is used to check how many genes in each cluster that have array
     * expression values.
     * @throws Exception
     */
    @Test
    public void checkClusterGenesOverlapping() throws Exception {
        String[] fileNames = new String[] {
                dirName + "NejmLogRatioZScore.txt",
                dirName + "GSE2034Filtered.txt"
        };
        Map<Integer, Set<String>> clusterToGenes = getClusterToGenes();
        String line = null;
        for (String fileName : fileNames) {
//            fu.setInput(fileName);
//            Set<String> genesInArray = new HashSet<String>();
//            line = fu.readLine();
//            while ((line = fu.readLine()) != null) {
//                String[] tokens = line.split("\t");
//                genesInArray.add(tokens[0]);
//            }
//            fu.close();
            Set<String> genesInArray = null;
            if (fileName.contains("Nejm"))
                genesInArray = getNejmGenes();
            else
                genesInArray = getGPL96Genes();
            System.out.println(fileName);
            int count = 0;
            for (Iterator<Integer> it = clusterToGenes.keySet().iterator(); it.hasNext();) {
                Integer cluster = it.next();
                Set<String> genes = clusterToGenes.get(cluster);
                Set<String> copy = new HashSet<String>(genes);
                copy.retainAll(genesInArray);
                double ratio = (double)copy.size() / genes.size();
                if (ratio >= 0.5)
                    count ++;
                System.out.println(cluster + "\t" + genes.size() + "\t" + 
                                   copy.size() + "\t" + ratio);
            }
            System.out.println("Clustering having 50% coverage: " + count);
        }   
    }
    
    
    /**
     * This method is used to check the overlapping between genes in the clusters and 
     * genes in the DNA array data sets.
     * @throws Exception
     */
    @Test
    public void checkOverlapping() throws Exception {
        String[] fileNames = new String[] {
                dirName + "NejmLogRatioZScore.txt",
                dirName + "GSE2034Filtered.txt"
        };
        Map<Integer, Set<String>> clusterToGenes = getClusterToGenes();
        Set<String> genesInClusters = new HashSet<String>();
        for (Set<String> set : clusterToGenes.values())
            genesInClusters.addAll(set);
        String line = null;
        System.out.println("Total genes in clusters: " + genesInClusters.size());
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            Set<String> genesInArray = new HashSet<String>();
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                genesInArray.add(tokens[0]);
            }
            fu.close();
            // How many genes in the array data set
            Set<String> copy = new HashSet<String>(genesInClusters);
            copy.removeAll(genesInArray);
            System.out.println(fileName);
            System.out.println("Total genes in array: " + genesInArray.size());
            System.out.println("Covered genes: " + (genesInArray.size() - copy.size()) / (double)genesInArray.size());
//            int c = 0;
//            for (String id : copy) {
//                System.out.println("\t" + id);
//                if (++c > 10)
//                    break;
//            }
        }
    }
    
    public Map<Integer, Set<String>> getClusterToGenes() throws Exception {
        Map<Integer, Set<String>> clusterToProteins = getClusterToProteins();
        // Check total covered proteins
        Set<String> totalProteins = new HashSet<String>();
        for (Set<String> set : clusterToProteins.values())
            totalProteins.addAll(set);
        System.out.println("Total proteins: " + totalProteins.size());
        Map<String, Set<String>> idToNames = new HibernateFIReader().generateProteinNameToAccession();
        Map<Integer, Set<String>> clusterToGenes = new HashMap<Integer, Set<String>>();
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        //Map<String, String> idToNames = uniAnalyzer.loadUniProtToGeneNameMap();
        for (Iterator<Integer> it = clusterToProteins.keySet().iterator(); it.hasNext();) {
            Integer index = it.next();
            Set<String> proteins = clusterToProteins.get(index);
            Set<String> genes = new HashSet<String>();
            for (String protein : proteins) {
                Set<String> names = idToNames.get(protein);
                if (names != null && names.size() == 1) {
                    // Use single mapped genes only
                    genes.add(names.iterator().next());
                }
                //String name = idToNames.get(protein);
                //genes.add(protein);
            }
            clusterToGenes.put(index, genes);
        }
        return clusterToGenes;
    }
    
    /**
     * This method is used to get proteins in clusters.
     * @return
     * @throws IOException
     */
    private Map<Integer, Set<String>> getClusterToProteins() throws IOException {
        // Load clusters
        String clusterFile = R3Constants.CLUSTER_RESULT_FILE;
        fu.setInput(clusterFile);
        Map<Integer, Set<String>> clusterToProteins = new HashMap<Integer, Set<String>>();
        int lastIndex = 119;
        int currentIndex = 1;
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Set<String> proteins = new HashSet<String>();
            for (String token : tokens)
                proteins.add(token);
            clusterToProteins.put(currentIndex,
                                  proteins);
            currentIndex ++;
            if (currentIndex > lastIndex)
                break;
        }
        fu.close();
        return clusterToProteins;
    }
    
    /**
     * This method is used to run z-normalization on a gene-wise based.
     * @throws IOException
     */
    @Test
    public void zNormalizeValues() throws IOException {
        String[] fileNames = new String[] {
                //dirName + "NejmLogRatio.txt",
//                dirName + "GSE2034FilteredNorm.txt",
                dirName + "NejmLogRatioNorm.txt"
//                dirName + "GSE4922FilteredOnSamples.txt"
        };
        for (String fileName : fileNames) {
            //zNormalizeValues(fileName);
            globalzNormalizeValues(fileName);
        }
    }
    
    protected Map<String, String> loadNejmPatientInfo() throws IOException {
        Map<String, String> sampleToType = new HashMap<String, String>();
        String fileName =  dirName + "Nejm_ClinicalData_Table.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            sampleToType.put("Sample " + tokens[1],
                             tokens[4]);
        }
        fu.close();
        return sampleToType;
    }
    
    protected Map<String, String> loadGSE4922PatientInfo() throws IOException {
        Map<String, String> sampleToType = new HashMap<String, String>();
        String fileName = dirName + "/GSE4922/Patient_info.txt";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            if (line.startsWith("GSM ID"))
                continue;
            String[] tokens = line.split("\t");
            if (tokens.length < 9)
                continue; // No patient information
            if (tokens[8].length() == 0)
                continue;
            String sample = tokens[0];
            int index = sample.indexOf("/");
            sample = sample.substring(0, index);
            sampleToType.put(sample, tokens[8]);
        }
        fu.close();
        return sampleToType;
    }
    
    private void zNormalizeValues(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        FileUtility outFu = new FileUtility();
        int index = fileName.lastIndexOf(".");
        String outFileName = fileName.substring(0, index) + "ZScore.txt";
        outFu.setOutput(outFileName);
        outFu.printLine(line);
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            builder.append(tokens[0]);
            List<String> zvalues = zNormalize(tokens);
            for (String value : zvalues)
                builder.append("\t").append(value);
            outFu.printLine(builder.toString());
            builder.setLength(0);
        }
        outFu.close();
    }
    
    private void globalzNormalizeValues(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get the overall average and sd from all values
        DescriptiveStatistics stat = new DescriptiveStatistics();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].length() == 0)
                    continue;
                stat.addValue(Double.parseDouble(tokens[i]));
            }
        }
        double mean = stat.getMean();
        double sd = stat.getStandardDeviation();
        fu.close();
        // Re-read
        fu.setInput(fileName);
        line = fu.readLine();
        FileUtility outFu = new FileUtility();
        int index = fileName.lastIndexOf(".");       
//        String outFileName = fileName.substring(0, index) + "ZScore_091409.txt";
        String outFileName = fileName.substring(0, index) + "GlobalZScore_070111.txt";
        outFu.setOutput(outFileName);
        outFu.printLine(line);
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            builder.append(tokens[0]);
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].length() == 0)
                    builder.append("\t");
                else {
                    double zscore = (Double.parseDouble(tokens[i]) - mean) / sd;
                    builder.append("\t").append(zscore);
                }
            }
            outFu.printLine(builder.toString());
            builder.setLength(0);
        }
        outFu.close();
    }
    
    private List<String> zNormalize(String[] tokens) {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        List<String> normalizedValues = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++) {
            String valueStr = tokens[i];
            if (valueStr.length() == 0)
                continue;
            double value = Double.parseDouble(valueStr);
            stat.addValue(value);
        }
        double mean = stat.getMean();
        double sd = stat.getStandardDeviation();
        for (int i = 1; i < tokens.length; i++) {
            String valueStr = tokens[i];
            if (valueStr.length() == 0) {
                normalizedValues.add("");
                continue;
            }
            double value = Double.parseDouble(valueStr);
            double newValue = (value - mean) / sd;
            normalizedValues.add(newValue + "");
        }
        return normalizedValues;
    }
    
    /**
     * This method is used to consolidate two files for GSE2034 and extract single mapped data points
     * to be analyzed.
     * @throws IOException
     */
    @Test
    public void consolidateGSE2034datasets() throws IOException {
        Map<String, List<String>> probesetToGenes = getGPL96ProbesetToGenes();
        System.out.println("Mapping to one gene probe sets: " + probesetToGenes.size());
//        String[] fileNames = new String[] {
//                dirName + "GSE2034_series_matrix-1.txt",
//                dirName + "GSE2034_series_matrix-2.txt"
//        };
        String[] fileNames = new String[] {
                dirName + "GSE4922/GSE4922-GPL96_series_matrix-1.txt",
                dirName + "GSE4922/GSE4922-GPL96_series_matrix-2.txt"
        };
        String line = null;
        List<String> sampleNames = new ArrayList<String>();
        Map<String, List<String>> probesetToValues = new HashMap<String, List<String>>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            while ((line = fu.readLine()) != null) {
                if (line.length() == 0)
                    continue;
                if (line.startsWith("\"ID_REF\"")) {
                    extractGSE2034Samples(line, sampleNames);
                    continue; // Title line
                }
                if (line.startsWith("!"))
                    continue; // comments
                String[] tokens = line.split("\t");
                String probeset = removeQuotationMarks(tokens[0]);
                if (!probesetToGenes.containsKey(probeset)) 
                    continue;
                List<String> values = probesetToValues.get(probeset);
                if (values == null) {
                    values = new ArrayList<String>();
                    probesetToValues.put(probeset, values);
                }
                extractGSE2034Values(tokens, values);
            }
            fu.close();
        }
        // Do a simple average for genes that can be mapped to multiple probesets
        Map<String, List<String>> geneToValues = new HashMap<String, List<String>>();
        // Reverse map from probeset to genes to genes to probeset
        Map<String, List<String>> geneToProbesets = reverseMap(probesetToGenes);
        for (Iterator<String> it = geneToProbesets.keySet().iterator(); it.hasNext();) {
            String gene = it.next();
            List<String> probesetList = geneToProbesets.get(gene);
            List<String> geneValues = averageValues(probesetList,
                                                    probesetToValues);
            geneToValues.put(gene, geneValues);
        }
        //String outFileName = dirName + "GSE2034Filtered.txt";
        String outFileName = dirName + "GSE4922Filtered.txt";
        outputConsolidateValues(geneToValues, 
                                sampleNames, 
                                "gene",
                                outFileName);
    }
    
    /**
     * This method is used to average several row values
     * @param probesetList
     * @param probesetToValues
     * @return
     */
    private List<String> averageValues(Collection<String> probesetList,
                                       Map<String, List<String>> probesetToValues) {
        if (probesetList.size() == 1) {
            String probeset = probesetList.iterator().next();
            List<String> values = probesetToValues.get(probeset);
            return new ArrayList<String>(values);
        }
        List<List<Double>> allData = new ArrayList<List<Double>>();
        // Peek
        List<String> values = null;
        for (String probeset : probesetList) {
            values = probesetToValues.get(probeset);
            if (values != null)
                break;
        }
        // Initialize the above array
        for (int i = 0; i < values.size(); i++)
            allData.add(new ArrayList<Double>());
        for (String probeset : probesetList) {
            values = probesetToValues.get(probeset);
            if (values == null)
                continue;
                //throw new IllegalStateException(probeset + " has no data values!");
            for (int i = 0; i < values.size(); i++) {
                List<Double> list = allData.get(i);
                String value = values.get(i);
                if (value.length() == 0)
                    list.add(null);
                else
                    list.add(new Double(value));
            }
        }
        List<String> average = new ArrayList<String>(allData.size());
        for (List<Double> colValues : allData) {
            double total = 0.0d;
            int c = 0;
            for (Double value : colValues) {
                if (value == null)
                    continue;
                c ++;
                total += value;
            }
            if (c == 0)
                average.add("");
            else
                average.add(total / c + "");
        }
        return average;
    }
    
    private Map<String, List<String>> reverseMap(Map<String, List<String>> map) {
        Map<String, List<String>> reverseMap = new HashMap<String, List<String>>();
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            List<String> values = map.get(key);
            for (String value : values) {
                List<String> values1 = reverseMap.get(value);
                if (values1 == null) {
                    values1 = new ArrayList<String>();
                    reverseMap.put(value, values1);
                }
                values1.add(key);
            }
        }
        return reverseMap;
    }
    
    /**
     * This method is used to remove quotation marks in the probeset ids.
     */
    private String removeQuotationMarks(String id) {
        int index1 = id.indexOf("\"");
        int index2 = id.lastIndexOf("\"");
        if (index1 == index2) {
            return id.substring(index1 + 1);
        }
        else
            return id.substring(index1 + 1, index2);
    }
    
    private void extractGSE2034Values(String[] tokens,
                                     List<String> values) {
        for (int i = 1; i < tokens.length; i++)
            values.add(tokens[i]);
    }
    
    protected void extractGSE2034Samples(String line,
                                       List<String> samples) {
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++)
            samples.add(removeQuotationMarks(tokens[i]));
    }
    
    /**
     * This helper method is used to get the map from probeset to genes for platform GPL96 in GEO.
     * @return
     * @throws IOException
     */
    private Map<String, List<String>> getGPL96ProbesetToGenes() throws IOException {
        String fileName = dirName + "GPL96.txt";
        fu.setInput(fileName);
        String line = null;
        Map<String, List<String>> probesetToGenes = new HashMap<String, List<String>>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            if (line.startsWith("ID"))
                continue;
            String[] tokens = line.split("\t");
            if (tokens.length < 11)
                continue;
            String gene = tokens[10];
            if (gene.length() == 0)
                continue;
            if (gene.contains("///"))
                continue; // A probeset can be mapped to more than one gene
            List<String> genes = probesetToGenes.get(tokens[0]);
            if (genes == null) {
                genes = new ArrayList<String>();
                probesetToGenes.put(tokens[0],
                                    genes);
            }
            genes.add(gene);
        }
        fu.close();
        return probesetToGenes;
    }
    
    /**
     * Three things are done here:
     * 1) Remove genes don't have enough data sets. 
     * 2) Average genes having more than one rows.
     * 3) Mapping gene names back to standard names.
     * @throws IOException
     */
    @Test
    public void normalizeArrayDataset() throws IOException {
        String sourceFile = dirName + "NejmLogRatio.txt";
        String targetFile = dirName + "NejmLogRatioNorm_011111.txt";
//        String sourceFile = dirName + "GSE2034Filtered.txt";
//        String targetFile = dirName + "GSE2034FilteredNorm.txt";
//        String sourceFile = dirName + "GSE4922FilteredOnSamples.txt";
//        String targetFile = dirName + "GSE4922FilteredOnSamplesNorm.txt";
        fu.setInput(sourceFile);
        String line = fu.readLine();
        List<String> samples = new ArrayList<String>();
        extractGSE2034Samples(line, samples);
        Map<String, List<String>> geneToValues = new HashMap<String, List<String>>();
        // Test shown that there are no gene duplicated in the file: 11189 in total
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            List<String> values = new ArrayList<String>();
            for (int i = 1; i < tokens.length; i++) {
                String token = tokens[i].trim();
                //if (token.length() > 0 && !token.equals("NA"))
                    values.add(token);
            }
            // Remove all genes that have null expression
            if (values.size() < samples.size())
                continue;
            geneToValues.put(gene, values);
        }
        System.out.println("Total genes: " + geneToValues.size());
        EntrezGeneAnalyzer analyzer = new EntrezGeneAnalyzer();
        Map<String, String> synonymToName = analyzer.getGeneSynonymToNameMap();
        Map<String, Set<String>> nameToSynonyms = new HashMap<String, Set<String>>();
        for (String gene : geneToValues.keySet()) {
            String name = synonymToName.get(gene);
            if (name == null)
                name = gene;
            Set<String> set = nameToSynonyms.get(name);
            if (set == null) {
                set = new HashSet<String>();
                nameToSynonyms.put(name, set);
            }
            set.add(gene);
        }
        Map<String, List<String>> nameToValues = new HashMap<String, List<String>>();
        for (String name : nameToSynonyms.keySet()) {
            Set<String> synonyms = nameToSynonyms.get(name);
            List<String> average = averageValues(synonyms, 
                                                 geneToValues);
            nameToValues.put(name, average);
        }
        System.out.println("After normalization: " + nameToValues.size());
//        outputConsolidateValues(nameToValues, 
//                                samples, 
//                                "gene", 
//                                targetFile);
    }
    
    /**
     * Use this method to generate one file for all log values for van de Vijver data set.
     * Only rows that have been validated and can be mapped to genes have been consolidated
     * together.
     * @throws IOException
     */
    @Test
    public void consolidateLogValueForNejmDataset() throws IOException {
        List<File> files = getNejmDataFiles();
        String outFileName = dirName + "NejmLogRatio.txt";
        // Want to get all sample names first
        List<String> sampleNames = new ArrayList<String>();
        String line = null;
        for (File file : files) {
            fu.setInput(file.getAbsolutePath());
            line = fu.readLine();
            String[] tokens = line.split("\t");
            for (String token : tokens) {
                if (token.length() > 0)
                    sampleNames.add(token);
            }
            fu.close();
        }
        System.out.println("Total samples: " + sampleNames.size());
        // Pick data from files
        // Want to get single genes only
        Set<String> geneSet = getSingleRowGenes(files.get(0));
        System.out.println("Total single row genes: " + geneSet.size());
        Map<String, List<String>> geneToValues = new HashMap<String, List<String>>();
        // Read logratio values from the original data set files.
        for (File file : files) {
            fu.setInput(file.getAbsolutePath());
            line = fu.readLine();
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                // Check if this gene should be in the output
                String geneName = tokens[1];
                if (!geneSet.contains(geneName)) 
                    continue;
                // Try to get value
                String logRatio = null;
                String flag;
                List<String> values = geneToValues.get(geneName);
                if (values == null) {
                    values = new ArrayList<String>();
                    geneToValues.put(geneName, values);
                }
                for (int i = 0; i < tokens.length; i++) {
                    if (i < 2)
                        continue; // escape the first two columns
                    int modulo = i % 5;
                    if (modulo == 2) // value
                        logRatio = tokens[i];
                    else if (modulo == 1) {
                        flag = tokens[i];
                        if (flag.equals("0")) {
                            values.add("");
                        }
                        else {
                            values.add(logRatio);
                        }
                    }
                }
            }
            fu.close();
        }
        outputConsolidateValues(geneToValues, 
                                sampleNames, 
                                "gene",
                                outFileName);
    }
    
    private void outputConsolidateValues(Map<String, List<String>> geneToValues,
                                         List<String> sampleNames,
                                         String firstColName,
                                         String outFileName) throws IOException {
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Print out the header line
        builder.append(firstColName);
        for (String sample : sampleNames)
            builder.append("\t").append(sample);
        outFu.printLine(builder.toString());
        // Output
        for (Iterator<String> it = geneToValues.keySet().iterator(); it.hasNext();) {
            String gene = it.next();
            List<String> values = geneToValues.get(gene);
            // Do a sanity check: if all values are empty, don't output it.
            boolean escape = true;
            for (String value : values) {
                if (!value.equals("")) {
                    escape = false;
                    break;
                }
            }
            if (escape)
                continue;
            builder.setLength(0);
            builder.append(gene);
            for (String value : values)
                builder.append("\t").append(value);
            outFu.printLine(builder.toString());
        }
        outFu.close();
    }
    
    /**
     * This help method is used to get the list of all original data files for van de Vijver
     * data set.
     * @return
     * @throws IOException
     */
    private List<File> getNejmDataFiles() throws IOException {
        String nejmDirName = this.dirName + "zipFiles295Samples";
        File dir = new File(nejmDirName);
        File[] files = dir.listFiles();
        List<File> fileList = new ArrayList<File>();
        for (File file : files) {
            if (file.getName().startsWith("Table"))
                fileList.add(file);
        }
        return fileList;
    }
    
    /**
     * This method is used to count how many genes in the van de Vijver data set.
     * @throws IOException
     */
    @Test
    public void analyzeNejmDataset() throws IOException {
        List<File> files = getNejmDataFiles();
        // Want to compare if any genes are duplicated displayed
        List<String> geneList = new ArrayList<String>();
        for (File file : files) {
            String fileName = file.getName();
            if (!fileName.startsWith("Table"))
                continue;
            grepNejmGenes(file, geneList);
            // Check for the first file only since the platform should be the same
            break;
        }
        Set<String> geneSet = new HashSet<String>(geneList);
        System.out.println("Genes in a list: " + geneList.size());
        System.out.println("Genes in a set: " + geneSet.size());
        Map<String, Integer> geneToNumber = InteractionUtilities.countTermUsageInList(geneList);
        // Want to find genes that have been listed once
        List<String> singleMappedGenes = new ArrayList<String>();
        for (Iterator<String> it = geneToNumber.keySet().iterator(); it.hasNext();) {
            String gene = it.next();
            Integer number = geneToNumber.get(gene);
            if (number == 1)
                singleMappedGenes.add(gene);
        }
        System.out.println("Single mapped genes: " + singleMappedGenes.size());
        //for (String gene : singleMappedGenes)
        //    System.out.println(gene);
    }
    
    @Test
    public void analyzeGSE4922Dataset() throws IOException {
        String[] fileNames = new String[] {
                dirName + "GSE4922/GSE4922-GPL96_series_matrix-1.txt",
                dirName + "GSE4922/GSE4922-GPL96_series_matrix-2.txt"
        };
        List<String> samples = new ArrayList<String>();
        List<String> dfsInfo = new ArrayList<String>();
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            String line = null;
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("!Sample_characteristics_ch1") &&
                    line.contains("DFS EVENT (0=censored;")) {
                    String[] tokens = line.split("\t");
                    for (int i = 1; i < tokens.length; i++) {
                        if (tokens[i].length() > 0)
                            dfsInfo.add(tokens[i]);
                    }
                }
                else if (line.startsWith("\"ID_REF\"")) {
                    String[] tokens = line.split("\t");
                    for (int i = 1; i < tokens.length; i++) 
                        samples.add(tokens[i]);
                }
            }
            fu.close();
        }
        System.out.println("Total samples: " + samples.size());
        System.out.println("Total samples having DFS info:" + dfsInfo.size());
        for (int i = 0; i < samples.size(); i++) {
            if (i > dfsInfo.size() - 1)
                System.out.println(samples.get(i));
            else
                System.out.println(samples.get(i) + ": " + dfsInfo.get(i));
        }
    }
    
    public Set<String> getNejmGenes() throws IOException {
        List<File> files = getNejmDataFiles();
        // Want to compare if any genes are duplicated displayed
        List<String> nejmGeneList = new ArrayList<String>();
        for (File file : files) {
            String fileName = file.getName();
            if (!fileName.startsWith("Table"))
                continue;
            grepNejmGenes(file, nejmGeneList);
            // Check for the first file only since the platform should be the same
            break;
        }
        return new HashSet<String>(nejmGeneList);
    }
    
    private void grepNejmGenes(File file,
                               List<String> geneList) throws IOException {
        fu.setInput(file.getAbsolutePath());
        // Top two lines are for annotation
        String line = fu.readLine();
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String geneName = tokens[1];
            if (geneName.length() > 1)
                geneList.add(geneName.toUpperCase());
        }
        fu.close();
    }
    
    private Set<String> getSingleRowGenes(File file) throws IOException {
        List<String> geneList = new ArrayList<String>();
        grepNejmGenes(file, geneList);
        Map<String, Integer> geneToNumber = InteractionUtilities.countTermUsageInList(geneList);
        Set<String> singleGenes = new HashSet<String>();
        for (Iterator<String> it = geneToNumber.keySet().iterator(); it.hasNext();) {
            String gene = it.next();
            Integer count = geneToNumber.get(gene);
            if (count == 1)
                singleGenes.add(gene);
        }
        return singleGenes;
    }
}
