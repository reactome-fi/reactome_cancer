/*
 * Created on May 21, 2015
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.cancer.driver.FICancerDriverPredictor;
import org.reactome.fi.pgm.FIMAFFileLoader;
import org.reactome.fi.pgm.FIPGMConfiguration;
import org.reactome.fi.pgm.FIPGMRunner;
import org.reactome.r3.DbNSFPAnalyzer;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;
import org.reactome.r3.util.UniProtProteinLengthHelper;

/**
 * This class is used to perform ICGC PanCancer Analysis.
 * @author gwu
 *
 */
public class ICGCPanCancerAnalyzer {
    private final String ICGC_PANCANCER_RESULT = "results/ICGC_PanCancer/2016_04/";
    //    private final String MAF_FILE_NAME = "datasets/ICGC/Santa_Cruz_Pilot/AnnotatedMAFs/Sanger.PCAWG_train2.annotated.snv_mnv.maf";
    //    private final String MAF_FILTERED_FILE_NAME = "datasets/ICGC/Santa_Cruz_Pilot/AnnotatedMAFs/Sanger.PCAWG_train2.annotated.snv_mnv.filter.genes.maf";
    private final String DATA_DIR = "datasets/ICGC/2016_04/";
    private final String MAF_FILE_NAME = DATA_DIR + "Barcelona_consensus.maf";
    private final String MAF_FILTERED_FILE_NAME = DATA_DIR + "Barcelona_consensus.filter.genes.maf";
    private final FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public ICGCPanCancerAnalyzer() {
    }
    
    @Test
    public void analyzeConnectionDegrees() throws Exception {
        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
        Set<String> fis = config.getFIs();
        final Map<String, Set<String>> proteinToPartners = InteractionUtilities.generateProteinToPartners(fis);
        List<String> proteins = new ArrayList<String>(proteinToPartners.keySet());
        Collections.sort(proteins, new Comparator<String>() {
            public int compare(String protein1, String protein2) {
                Set<String> partners1 = proteinToPartners.get(protein1);
                Set<String> partners2 = proteinToPartners.get(protein2);
                return partners2.size() - partners1.size();
            }
        });
        for (int i = 0; i < proteins.size(); i++) {
            String protein = proteins.get(i);
            System.out.println(protein + "\t" + proteinToPartners.get(protein).size());
            if (i == 100)
                break;
        }
    }
    
    @Test
    public void convertTextToARFF() throws Exception {
        String src = DATA_DIR + "PGM_FI_Inference_Results_071416_31_19_Annot.txt";
        String target = DATA_DIR + "PGM_FI_Inference_Results_071416_31_19_Annot.arff";
        fu.setInput(src);
        fu.setOutput(target);
        fu.printLine("@RELATION icgc");
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        for (String token : tokens)
            fu.printLine("@ATTRIBUTE " + token + " NUMERIC");
        fu.printLine("@DATA");
        while ((line = fu.readLine()) != null) {
            String line1 = line.replaceAll("\t", ",");
            fu.printLine(line1);
        }
        fu.close();
    }

    @Test
    public void convertToMatrix() throws Exception {
        FIPGMRunner runner = new FIPGMRunner();
//        String fileName = DATA_DIR + "PGM_FI_Inference_Results_060216.txt";
        String fileName = DATA_DIR + "PGM_FI_Inference_Results_061316.txt";
        fileName = DATA_DIR + "PGM_FI_Inference_Results_071216.txt";
        fileName = DATA_DIR + "PGM_FI_Inference_Results_071316_1.txt";
        
        int index = fileName.lastIndexOf(".");
        String tmpFileName = fileName.substring(0, index);
        
        String targetFileName = tmpFileName + "_LogRatio.txt";
        runner.convertToMatrix(fileName,
                               targetFileName,
                               5);
        
        targetFileName = tmpFileName + "_Observation.txt";
        runner.convertToMatrix(fileName,
                               targetFileName,
                               7);
        
        // Generate a two_sum file name
        Map<String, Double> geneToLogRatio = openGeneToScore(tmpFileName + "_LogRatio.txt");
        Map<String, Double> geneToObservation = openGeneToScore(tmpFileName + "_Observation.txt");
        targetFileName = tmpFileName + "Two_Sums.txt";
        fu.setOutput(targetFileName);
        fu.printLine("Gene\tLogRatio\tObservation");
        for (String gene : geneToLogRatio.keySet()) {
            fu.printLine(gene + "\t" +
                         geneToLogRatio.get(gene) + "\t" + 
                         geneToObservation.get(gene));
        }
        fu.close();
        annotateResultFiles(targetFileName);
    }
    
    private Map<String, Double> openGeneToScore(String fileName) throws IOException {
        fu.setInput(fileName);
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToScore.put(tokens[0], new Double(tokens[1]));
        }
        fu.close();
        return geneToScore;
    }
    
    /**
     * Perform MRF analysis.
     * @throws Exception
     */
    @Test
    public void performMRFAnalysis() throws Exception {
        FIPGMRunner runner = new FIPGMRunner();
        runner.runInferenceBasedOnDataTypes();
    }

    /**
     * Check mutation scores in one single sample after loading a maf file.
     * @throws IOException
     */
    @Test
    public void checkMutationScoresInOneSample() throws IOException {
        String fileName = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.maf";
        MATFileLoader matFileLoader = new MATFileLoader();
        // Set as null to take the whole texts in the sample columns
        matFileLoader.setSampleNameLength(null);
        Map<String, Map<String, Float>> sampleToGeneToScore = matFileLoader.loadSampleToGeneToFIScore(fileName,
                                                                                                     "MetaLR_RankScore");
        System.out.println("Total samples: " + sampleToGeneToScore.size());
        System.out.println("Sample\tMutated_Genes");
        for (String sample : sampleToGeneToScore.keySet()) {
            Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
            System.out.println(sample + "\t" + geneToScore.size());
        }
        // This sample has an extremely large size of mutated genes
        String sample = "bcf858fd-cc3b-4fde-ab10-eb96216f4366";
        Map<String, Float> geneToScore = sampleToGeneToScore.get(sample);
        System.out.println("\n" + sample);
        System.out.println("Gene\tScore");
        for (String gene : geneToScore.keySet()) {
            Float score = geneToScore.get(gene);
            System.out.println(gene + "\t" + score);
        }
    }
    
    @Test
    public void checkdbNSFPScores() throws IOException {
        String fileName = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.maf";
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int size = tokens.length;
        // We want to collect values in the last column
        DescriptiveStatistics stat = new DescriptiveStatistics();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (tokens.length < size)
                continue;
            if (tokens[size - 1].length() == 0)
                continue;
            stat.addValue(new Double(tokens[size - 1]));
        }
        fu.close();
        System.out.println("Maximum: " + stat.getMax());
        System.out.println("Minium: " + stat.getMin());
        System.out.println("Mean: " + stat.getMean());
        System.out.println("Variance: " + stat.getVariance());
        System.out.println("Medium: " + stat.getPercentile(50.0d));
        System.out.println("Total data points: " + stat.getN());
        // Want to print out all scores
        String output = DATA_DIR + "AnnotatedScores.txt";
        fu.setOutput(output);
        for (int i = 0; i < stat.getN(); i++) {
            double value = stat.getElement(i);
            fu.printLine(value + "");
        }
        fu.close();
    }
    
    @Test
    public void runAnnotateResultFiles() throws Exception {
//        annotateResultFiles(DATA_DIR + "PGM_FI_Inference_Results_071416.txt");
        annotateResultFiles(DATA_DIR + "PGM_FI_Inference_Results_071816_TwoScores.txt");
    }
    
    private void annotateResultFiles(String inFileName) throws Exception {
        Map<String, Integer> geneToLength = new UniProtProteinLengthHelper().loadGeneToProteinLength();
        FIPGMConfiguration config = FIPGMConfiguration.getConfig();
        Set<String> fis = config.getFIs();
        Map<String, Set<String>> proteinToPartners = InteractionUtilities.generateProteinToPartners(fis);
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> cgcGenes = helper.loadCancerCensusGenes();
        Set<String> allCGCGenes = helper.loadCancerCensusGenes(null, null);
        
        int index = inFileName.lastIndexOf(".");
        String outFileName = inFileName.substring(0, index) + "_Annot.txt";
        fu.setInput(inFileName);
        fu.setOutput(outFileName);
        
        String line = fu.readLine();
        fu.printLine(line + "\tLength\tDegree\tCGC\tAll_CGC\tLogRatio/Length\tObservation/Length");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Integer length = geneToLength.get(tokens[0]);
            if (length == null)
                continue; // Just ignore them!
            Set<String> partners = proteinToPartners.get(tokens[0]);
            if (partners == null) {
                System.err.println(tokens[0] + " has no partners!");
                continue;
            }
            Double logRatioToLength = null;
            if (length != null)
                logRatioToLength = new Double(tokens[1]) / length;
            Double observationToLength = null;
            if (length != null)
                observationToLength = new Double(tokens[2]) / length;
            fu.printLine(line + "\t" + length + "\t" + partners.size() + "\t" + 
                         (cgcGenes.contains(tokens[0]) ? 1 : 0) + "\t" + 
                         (allCGCGenes.contains(tokens[0]) ? 1 : 0) + "\t" + 
                         logRatioToLength + "\t" + 
                         observationToLength);
        }
        
        fu.close();
    }
    
    /**
     * Divide the mutation scores by protein lengths to normalize them.
     * @throws IOException
     */
    @Test
    public void normalizeMutationScores() throws IOException {
        String srcFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.maf";
//        String targetFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.normalized.maf";
        String targetFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.normalized.simple.maf";
        fu.setInput(srcFile);
        fu.setOutput(targetFile);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        // Score is the last column
        Map<String, Integer> proteinToLength = new UniProtProteinLengthHelper().loadGeneToProteinLength();
        fu.printLine(line + "\tNormalizedScore\tLog_NormalizedScore");
        int annotatedLines = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < headers.length) {
//                fu.printLine(line);
                continue; // No score
            }
            Integer length = proteinToLength.get(tokens[0]);
            if (length == null) {
//                fu.printLine(line);
                continue;
            }
            Double score = new Double(tokens[headers.length - 1]);
            Double nScore = score / length;
            Double logNScore = Math.log(nScore);
            fu.printLine(line + "\t" + nScore + "\t" + logNScore);
            annotatedLines ++;
        }
        fu.close();
        System.out.println("Total annotated lines: " + annotatedLines);
    }
    
    @Test
    public void generateSumOfScoreFile() throws Exception {
        String mafFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.071816.maf";
        String scoreCol = "MetaLR_score";
        String mafSumFile = DATA_DIR + "ICGC_Score_Sum_071816.txt";
        
        FIMAFFileLoader mafFileLoader = new FIMAFFileLoader();
        Map<String, Map<String, Float>> sampleToGeneToScore = mafFileLoader.loadSampleToGeneToFIScore(mafFile, 
                                                                                                      scoreCol);
        // Get all genes
        Set<String> allGenes = new HashSet<String>();
        for (Map<String, Float> geneToScore : sampleToGeneToScore.values())
            allGenes.addAll(geneToScore.keySet());
        System.out.println("Total genes: " + allGenes.size());
        
        fu.setOutput(mafSumFile);
        fu.printLine("Gene\tScoreSum\tScoreSum/ProteinLength");
        Map<String, Integer> proteinToLength = new UniProtProteinLengthHelper().loadGeneToProteinLength();
        for (String gene : allGenes) {
            Integer length = proteinToLength.get(gene);
            if (length == null)
                continue; // Just escape it
            Double totalScore = 0.0d;
            for (Map<String, Float> geneToScore : sampleToGeneToScore.values()) {
                Float score = geneToScore.get(gene);
                if (score == null)
                    continue;
                totalScore += score;
            }
            fu.printLine(gene + "\t" + 
                         totalScore + "\t" + 
                         totalScore / length);
        }
        fu.close();
    }
    
    @Test
    public void annotateMAFFileWithdbNSFP() throws IOException {
        String scoreFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.output.txt";
        String outFileName = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.071819.maf";
        DbNSFPAnalyzer helper = new DbNSFPAnalyzer();
        helper.annotateMAFFileWithdbNSFP(MAF_FILE_NAME,
                                         outFileName,
                                         scoreFile);
    }
    
    @Test
    public void checkdbNSFPFiles() throws IOException {
        String srcFile = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.output.txt";
        fu.setInput(srcFile);
        String line = fu.readLine();
        fu.close();
        String[] tokens = line.split("\t");
        for (int i = 0; i < tokens.length; i++) {
            System.out.println(i + "\t" + tokens[i]);
        }
    }
    
    @Test
    public void generateInputFordbNSFP() throws IOException {
        DbNSFPAnalyzer analyzer = new DbNSFPAnalyzer();
        String outputFileName = DATA_DIR + "Barcelona_consensus.filter.genes.dbNSFP.input.txt";
        analyzer.generateInputFordbNSFPFromMAF(MAF_FILTERED_FILE_NAME,
                                               outputFileName);
    }
    
    @Test
    public void compareGeneLists() throws IOException {
        String[] fileNames = new String[] {
                "ICGCPanCancerScores_Network_Only_061115.txt",
                "ICGCPanCancerScores_NoMutation_052815.txt",
                "ICGCPanCancerScores_NoMutSig_061115.txt"
        };
        List<List<String>> lists = new ArrayList<List<String>>();
        double cutoff = 0.50d;
        for (String fileName : fileNames) {
            fu.setInput(ICGC_PANCANCER_RESULT + fileName);
            String line = fu.readLine();
            List<String> list = new ArrayList<String>();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                double score = new Double(tokens[tokens.length - 1]);
                if (score < cutoff)
                    continue;
                list.add(tokens[0]);
            }
            fu.close();
            lists.add(list);
        }
        for (int i = 0; i < lists.size() - 1; i++) {
            List<String> list1 = lists.get(i);
            System.out.println("List1 size: " + list1.size());
            for (int j = i + 1; j < lists.size(); j++) {
                List<String> list2 = lists.get(j);
                System.out.println("List2 size: " + list2.size());
                Set<String> shared = InteractionUtilities.getShared(list1, list2);
                System.out.println("Shared: " + shared.size());
                List<String> list2Copy = new ArrayList<String>(list2);
                list2Copy.removeAll(list1);
                System.out.println("Not shared in list2: " + list2Copy.size());
                for (String gene : list2Copy)
                    System.out.println(gene);
                List<String> list1Copy = new ArrayList<String>(list1);
                list1Copy.removeAll(shared);
                System.out.println("Not shared in list1: " + list1Copy.size());
                for (String gene : list1Copy)
                    System.out.println(gene);
            }
        }
    }
    
    /**
     * Two methods in this method are in class CancerDriverInstanceGenerator. Comment this
     * method out and should refactor later on.
     * @throws Exception
     */
//    @Test
//    public void calculateMutationExclusivity() throws Exception {
//        String[] genes = new String[]{"FGFR2", "FGFR4", "FGFR3", "KLB"};
////        genes = new String[]{"TMEM204", "FLT1", "FLT4", "GPR116", "KDR"};
////        genes = new String[]{"VHL", "HINT1", "RNF139", "OTUD7A"};
////        genes = new String[]{"HSD3B2", "CGA", "HSD3B1", "HSD17B3"};
//        String mafFileName = getMAFFileName();
//        MATFileLoader fileLoader = new MATFileLoader();
//        
//        Map<String, Set<String>> sampleToGenes = fileLoader.loadSampleToGenes(mafFileName,
//                                                                              "Protein_Change",
//                                                                              null,
//                                                                              new UniProtAnalyzer().loadGeneToProteinLength());
//        Set<String> hyperSamples = getHyperMutatedSamples();
//        System.out.println("Total samples in MAF: " + sampleToGenes.size());
//        sampleToGenes.keySet().removeAll(hyperSamples);
//        System.out.println("After removing hyper mutated samples: " + sampleToGenes.size());
//        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
//        // Get the total samples
//        Set<String> totalSamples = new HashSet<String>();
//        for (String gene : genes)
//            totalSamples.addAll(geneToSamples.get(gene));
//        FisherExact fisher = new FisherExact(sampleToGenes.size());
//        for (int i = 0; i < genes.length - 1; i++) {
//            String gene1 = genes[i];
//            Set<String> samples1 = geneToSamples.get(gene1);
//            for (int j = i + 1; j < genes.length; j++) {
//                String gene2 = genes[j];
//                Set<String> samples2 = geneToSamples.get(gene2);
//                Set<String> shared = InteractionUtilities.getShared(samples1, samples2);
//                Set<String> b = new HashSet<String>(samples2);
//                b.removeAll(shared);
//                Set<String> c = new HashSet<String>(samples1);
//                c.removeAll(shared);
//                int a = totalSamples.size() - shared.size() - b.size() - c.size();
//                double leftP = fisher.getLeftTailedP(a, b.size(), c.size(), shared.size());
//                double rightP = fisher.getRightTailedP(a, b.size(), c.size(), shared.size());
//                System.out.println(a + "\t" + b.size() + "\t" + c.size() + "\t" + shared.size());
//                System.out.println(gene1 + "\t" + gene2 + "\t" + leftP + "\t" + rightP);
//            }
//        }
//        List<String> samples = new ArrayList<String>(totalSamples);
//        Collections.sort(samples);
//        StringBuilder builder = new StringBuilder();
//        builder.append("Sampe");
//        for (String gene : genes)
//            builder.append("\t").append(gene);
//        System.out.println(builder.toString());
//        builder.setLength(0);
//        int singleMutationCount = 0;
//        for (String sample : samples) {
//            builder.append(sample);
//            int count = 0;
//            for (String gene : genes) {
//                builder.append("\t");
//                Set<String> geneSamples = geneToSamples.get(gene);
//                if (geneSamples.contains(sample)) {
//                    builder.append("&");
//                    count ++;
//                }
//            }
//            System.out.println(builder.toString());
//            builder.setLength(0);
//            if (count == 1)
//                singleMutationCount ++;
//        }
//        System.out.println("SingleMutationCount: " + singleMutationCount + " (" + (double) singleMutationCount / samples.size() + ")");
//        System.out.println("\tRandom test");
//        System.out.println("TotalSamples\tTotalSingleMutationSamples\tRatio");
//        for (int i = 0; i < 100; i++) {
//            List<Set<String>> randomSets = new ArrayList<Set<String>>();
//            Set<String> totalRandomSamples = new HashSet<String>();
//            for (String gene : genes) {
//                Set<String> geneSamples = geneToSamples.get(gene);
//                Set<String> randomSamples = MathUtilities.randomSampling(sampleToGenes.keySet(), 
//                                                                         geneSamples.size());
//                randomSets.add(randomSamples);
//                totalRandomSamples.addAll(randomSamples);
//            }
//            singleMutationCount = 0;
//            for (String sample : totalRandomSamples) {
//                int count = 0;
//                for (Set<String> randomSet : randomSets) {
//                    if (randomSet.contains(sample))
//                        count ++;
//                }
//                if (count == 1)
//                    singleMutationCount ++;
//            }
//            System.out.println(totalRandomSamples.size() + "\t" +
//                               singleMutationCount + "\t" + 
//                               (double)singleMutationCount / totalRandomSamples.size());
//        }
//    }
    
    @Test
    public void performMCLClusteringBasedLogistScore() throws Exception {
        String logisticResultFile = ICGC_PANCANCER_RESULT + "ICGCPanCancerScores_NoMutation_052815.txt";
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        fu.setInput(logisticResultFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToScore.put(tokens[0], new Double(tokens[tokens.length - 1]));
        }
        fu.close();
        
        MCLClusterWrapper clusterHelper = new MCLClusterWrapper();
        clusterHelper.setInflation(12.5d);
        List<Set<String>> clusters = clusterHelper.mclClusterForGeneScores(geneToScore,
                                                                           false,
                                                                           null,
                                                                           null);
        System.out.println("Total clusters: " + clusters.size());
        DescriptiveStatistics stat = new DescriptiveStatistics();
        System.out.println("ClusterIndex\tSize\tMedian\tMean\tGenes");
        for (int i = 0; i < clusters.size(); i++) {
            Set<String> cluster = clusters.get(i);
            if (cluster.size() < 3)
                continue;
            // Get the average score
            for (String gene : cluster) {
                Double score = geneToScore.get(gene);
                stat.addValue(score);
            }
            System.out.println(i + "\t" + 
                               cluster.size() + "\t" + 
                               stat.getPercentile(50.0d) + "\t" +
                               stat.getMean() + "\t" +
                               cluster);
            stat.clear();
        }
    }
    
    /**
     * Choosing one or more genes to link MutSig significant genes together.
     * @throws Exception
     */
    @Test
    public void increaseModularityByAddingGenes() throws Exception {
        //String fileName = ICGC_PANCANCER_RESULT + "ICGCPanCancerScores_052615.txt";
        String fileName = ICGC_PANCANCER_RESULT + "ICGCPanCancerScores_NoMutation_052815.txt";
        Set<String> seedGenes = new HashSet<String>();
        final Map<String, Double> geneToScore = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double score = new Double(tokens[8]);
            geneToScore.put(tokens[0], score);
            if (score >= 0.50d)
                seedGenes.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total seed genes: " + seedGenes.size());
        SpectralPartitionNetworkCluster clusterEngine = new SpectralPartitionNetworkCluster();
        Set<String> allFIs = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fisInSeedGenes = InteractionUtilities.getFIs(seedGenes, allFIs);
        List<Set<String>> clusters = clusterEngine.cluster(fisInSeedGenes);
        double modularity = clusterEngine.calculateModualarity(clusters,
                                                               fisInSeedGenes);
        System.out.println("Modularity: " + modularity);
        List<String> sortedGenes = new ArrayList<String>(geneToScore.keySet());
        sortedGenes.removeAll(seedGenes);
        Collections.sort(sortedGenes, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Double score1 = geneToScore.get(gene1);
                Double score2 = geneToScore.get(gene2);
                return score2.compareTo(score1);
            }
        });
        
        double current = modularity;
        Set<String> candidateGenes = new HashSet<String>(seedGenes);
        for (int i = 0; i < sortedGenes.size(); i++) {
            String gene = sortedGenes.get(i);
            Double score = geneToScore.get(gene);
//            System.out.println("Current score: " + score);
//            if (score < 0.07)
//                break;
            if (score < 0.125)
                break;
            candidateGenes.add(gene);
            fisInSeedGenes = InteractionUtilities.getFIs(candidateGenes, allFIs);
            clusters = clusterEngine.cluster(fisInSeedGenes);
            modularity = clusterEngine.calculateModualarity(clusters, fisInSeedGenes);
            if (modularity > current) {
                System.out.println(score + "\t" + 
                                   candidateGenes.size() + "\t" + 
                                   fisInSeedGenes.size() + "\t" +
                                   modularity);
//                System.out.println("Current Modularity: " + modularity);
//                System.out.println("Total candidate genes: " + candidateGenes.size());
//                System.out.println("Total FIs: " + fisInMutSigGenes.size());
                if (modularity > 0.70)
                    break;
                current = modularity;
                continue;
            }
            candidateGenes.remove(gene);
        }
        for (String gene : candidateGenes)
            System.out.println(gene);
    }
    
    /**
     * This method is used to check mutsigcv results for the ICGC pancer mutation file.
     * @throws Exception
     */
    @Test
    public void checkMutSigCVResults() throws Exception {
        String fileName = ICGC_PANCANCER_RESULT + "mutsigcv_results_052115_sig_genes.txt";
        double qvalueCutoff = 0.25d;
        Set<String> mutSigCVGenes = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            double qvalue = new Double(tokens[tokens.length - 1]);
            if (qvalue < qvalueCutoff)
                mutSigCVGenes.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total MutSigCV genes: " + mutSigCVGenes.size());
        // Check overlapping with known CGC genes
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> cgcGenes = helper.loadCancerCensusGenes();
        Set<String> shared = InteractionUtilities.getShared(mutSigCVGenes, cgcGenes);
        System.out.println("Shared with CGC genes: " + shared.size() + 
                           " (" + ((double)shared.size() / mutSigCVGenes.size()) + ")");
    }
    
    @Test
    public void checkSamplesInFilteredMAFFile() throws IOException {
        MATFileLoader fileLoader = new MATFileLoader();
        Map<String, Integer> geneToLength = new UniProtProteinLengthHelper().loadGeneToProteinLength();
//        Map<String, Set<String>> sampleToGenes = fileLoader.loadSampleToGenes(MAF_FILTERED_FILE_NAME,
//                                                                              "Protein_Change", 
//                                                                              null, 
//                                                                              geneToLength);
        // The Barcelona file has not aa change annotated
        Map<String, Set<String>> sampleToGenes = fileLoader.loadSampleToGenes(MAF_FILTERED_FILE_NAME,
                                                                              null, 
                                                                              null, 
                                                                              geneToLength);
        System.out.println("Total samples: " + sampleToGenes.size());
        // Some simple stats
        SummaryStatistics stats = new SummaryStatistics();
        String fileName = ICGC_PANCANCER_RESULT + "SampleToMutatedGeneNumber.txt";
        fu.setOutput(fileName);
        fu.printLine("Sample\tMutatedGenes");
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            stats.addValue(genes.size());
            fu.printLine(sample + "\t" + genes.size());
        }
        fu.close();
        System.out.println("Mean: " + stats.getMean() + ", sd: " + stats.getStandardDeviation());
//        List<String> sampleList = new ArrayList<String>(sampleToGenes.keySet());
//        Collections.sort(sampleList);
//        for (String sample : sampleList)
//            System.out.println(sample);
        Map<String, Set<String>> geneToSamples = InteractionUtilities.switchKeyValues(sampleToGenes);
        System.out.println("Total genes: " + geneToSamples.size());
        stats.clear();
        fileName = ICGC_PANCANCER_RESULT + "MutatedGeneToSampleNumber.txt";
        fu.setOutput(fileName);
        fu.printLine("MutatedGene\tSampleNumber");
        for (String gene : geneToSamples.keySet()) {
            Set<String> samples = geneToSamples.get(gene);
            stats.addValue(samples.size());
            fu.printLine(gene + "\t" + samples.size());
        }
        fu.close();
        System.out.println("Mean: " + stats.getMean() + ", sd: " + stats.getStandardDeviation());
    }
    
    /**
     * Remove lines with unknown and mutations inside introns.
     * @throws IOException
     */
    @Test
    public void filterMAFFileToGenesOnly() throws IOException {
        fu.setInput(MAF_FILE_NAME);
        String outFileName = MAF_FILTERED_FILE_NAME;
        fu.setOutput(outFileName);
        String line = fu.readLine();
        fu.printLine(line);
        int totalLines = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].equals("Unknown") || tokens[4].equals("Intron") || tokens[4].equals("RNA"))
                continue;
            if (tokens.length > 22) {
                String codonChange = tokens[20];
                String proteinChange = tokens[21];
                if (codonChange.length() == 0 && proteinChange.length() == 0)
                    continue;
            }
            fu.printLine(line);
            totalLines ++;
        }
        fu.close();
        System.out.println("Total lines in the new file: " + totalLines);
    }
    
    @Test
    public void checkChromosesInMAFFile() throws IOException {
        fu.setInput(MAF_FILTERED_FILE_NAME);
        String line = fu.readLine();
        List<String> headers = Arrays.asList(line.split("\t"));
        int chrIndex = headers.indexOf("Chromosome");
        System.out.println("Chromosome index: " + chrIndex);
        Set<String> chromosomes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            chromosomes.add(tokens[chrIndex]);
        }
        fu.close();
        List<String> list = new ArrayList<String>(chromosomes);
        Collections.sort(list);
        System.out.println("Total chromosomes: " + list.size());
        for (String chr : chromosomes)
            System.out.println(chr);
    }
    
    @Test
    public void checkMAFFile() throws IOException {
        fu.setInput(MAF_FILE_NAME);
        String line = fu.readLine();
        System.out.println(line);
        int count = 0;
        while ((line = fu.readLine()) != null) {
            System.out.println(line);
//            String[] tokens = line.split("\t");
//            if (tokens[0].equals("Unknown"))
//                continue;
//            String codonChange = tokens[20];
//            String proteinChange = tokens[21];
//            if (codonChange.length() == 0 && proteinChange.length() == 0) {
//                System.out.println(line);
//                continue;
//            }
//            System.out.println(tokens[0] + "\t" + codonChange + "\t" + proteinChange);
            count ++;
            if (count == 100)
                break;
        }
        fu.close();
    }
    
    @Test
    public void overlapAnalysis() throws IOException {
        String fileName = DATA_DIR + "PGM_FI_Inference_Results_071416_31_19_Annot.txt";
        fileName = DATA_DIR + "PGM_FI_Inference_Results_071416_26_24_Annot.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<GeneInfo> list = new ArrayList<ICGCPanCancerAnalyzer.GeneInfo>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            GeneInfo info = new GeneInfo();
            info.gene = tokens[0];
            info.logRatio = new Double(tokens[1]);
            info.normalizedScore = new Double(tokens[2]);
            info.degree = new Integer(tokens[4]);
            info.score = new Double(tokens[7]);
            list.add(info);
        }
        fu.close();
        
        Set<String> cgcGenesInTop = new HashSet<String>();
        FICancerDriverPredictor helper = new FICancerDriverPredictor();
        Set<String> cgcGenes = helper.loadCancerCensusGenes();
//        cgcGenes = helper.loadCancerCensusGenes(null, null);
        int top = 100;
        Set<String> allTopGenes = new HashSet<String>();
        
        Comparator<GeneInfo> sorter = new Comparator<ICGCPanCancerAnalyzer.GeneInfo>() {
            public int compare(GeneInfo info1, GeneInfo info2) {
                return info2.logRatio.compareTo(info1.logRatio);
            }
        };
        List<String> topGenes = getTopGenes(list, top, sorter);
        List<String> logRatioTop = new ArrayList<String>(topGenes);
        allTopGenes.addAll(topGenes);
        
        sorter = new Comparator<ICGCPanCancerAnalyzer.GeneInfo>() {
            public int compare(GeneInfo info1, GeneInfo info2) {
                return info2.normalizedScore.compareTo(info1.normalizedScore);
            }
        };
        topGenes = getTopGenes(list, top, sorter);
        List<String> normScoreTop = new ArrayList<String>(topGenes);
        allTopGenes.addAll(topGenes);
        System.out.println("LogRatioTop: " + logRatioTop.size());
        System.out.println("normScoreTop: " + normScoreTop.size());
        System.out.println("Shared: " + InteractionUtilities.getShared(normScoreTop, logRatioTop).size());
        Set<String> logRatioCGC = InteractionUtilities.getShared(logRatioTop, cgcGenes);
        Set<String> normScoreCGC = InteractionUtilities.getShared(normScoreTop, cgcGenes);
        System.out.println("logRatioCGC: " + logRatioCGC.size());
        System.out.println("normScoreCGC: " + normScoreCGC.size());
        System.out.println("Shared: " + InteractionUtilities.getShared(logRatioCGC, normScoreCGC).size());
        
        sorter = new Comparator<ICGCPanCancerAnalyzer.GeneInfo>() {
            public int compare(GeneInfo info1, GeneInfo info2) {
                return info2.degree.compareTo(info1.degree);
            }
        };
        topGenes = getTopGenes(list, top, sorter);
        allTopGenes.addAll(topGenes);
        
        sorter = new Comparator<ICGCPanCancerAnalyzer.GeneInfo>() {
            public int compare(GeneInfo info1, GeneInfo info2) {
                return info2.score.compareTo(info1.score);
            }
        };
        topGenes = getTopGenes(list, top, sorter);
        System.out.println("scoreTop: " + InteractionUtilities.getShared(cgcGenes, topGenes).size());
        allTopGenes.addAll(topGenes);
        
        System.out.println("Total top genes: " + allTopGenes.size());
        System.out.println("Total cgcGeneInTop: " + InteractionUtilities.getShared(allTopGenes, cgcGenes).size());
    }
    
    private List<String> getTopGenes(List<GeneInfo> list,
                                     int top,
                                     Comparator<GeneInfo> sorter) {
        Collections.sort(list, sorter);
        List<String> rtn = new ArrayList<String>();
        for (int i = 0; i < top; i++) {
            GeneInfo info = list.get(i);
            rtn.add(info.gene);
        }
        return rtn;
    }
    
    private class GeneInfo {
        
        String gene;
        Double logRatio;
        Double normalizedScore;
        Integer degree;
        Double score;
        
    }
    
}
