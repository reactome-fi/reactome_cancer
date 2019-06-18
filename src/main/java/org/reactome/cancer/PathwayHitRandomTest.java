/*
 * Created on Jun 6, 2007
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.junit.Test;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.genome.RefSeqInfo;
import org.reactome.genome.Segment;
import org.reactome.r3.TopicAnalyzer;
import org.reactome.r3.UCSCDataAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;


/**
 * This class is used to do random test for pathway hits for a list of genes.
 * TODO: This class needs to be revised
 * @author guanming
 *
 */
public class PathwayHitRandomTest {
//    private PathwayBasedAnnotator annotator;
//    private Map<String, Set<String>> idToTopics;
//    private Map<String, Double> topicToRatio;
//    private Map<String, Integer> topicToIdNumber;
//    
//    public PathwayHitRandomTest() throws IOException {
//        setUp();
//    }
//    
//    /**
//     * Initialize all necessary properties for doing test.
//     */
//    private void setUp() throws IOException {
//        annotator = new PathwayBasedAnnotator();
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        idToTopics = topicAnalyzer.generateIdToTopicWithHop(0);
//        int total = idToTopics.size();
//        topicToIdNumber = topicAnalyzer.countProteinsInTopics(idToTopics);
////        topicToRatio = annotator.calculateTopicToRatio(total, topicToIdNumber);
//    }
//    
//    @Test
//    public void singleRandomTest() throws IOException {
//        // Force to use a random seed
//        RandomData randomData = new RandomDataImpl(new JDKRandomGenerator());
//        Set<String> set = new HashSet<String>(idToTopics.keySet());
//        Object[] sampledIds = randomData.nextSample(set, 1120);
//        List<String> sample = new ArrayList<String>();
//        for (Object obj : sampledIds) {
//            sample.add(obj.toString());
//        }
//        annotator.topicAnnotate(sample, System.out);
//    }
//    
//    public void singleRandomTest(Collection<String> totalIds,
//                                 int geneNumberInSample,
//                                 PrintStream ps) throws IOException {
//        // Force to use a random seed
//        RandomData randomData = new RandomDataImpl(new JDKRandomGenerator());
//        Object[] sampledIds = randomData.nextSample(totalIds, geneNumberInSample);
//        List<String> sample = new ArrayList<String>();
//        for (Object obj : sampledIds) {
//            sample.add(obj.toString());
//        }
//        annotator.topicAnnotate(sample, System.out);
//    }
//    
//    @Test
//    public void randomTestForGeneNames() throws IOException {
//        setUp();
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        Set<String> geneNames = ucscAnalyzer.loadHumanGeneNames();
//        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
//        double threshold = 1.0E-6d;
//        List<String> hitTopics = new ArrayList<String>();
//        for (int i = 0; i < 1000; i++) {
//            // 3900 is the total gene names (nonredudant)
//            // 261 is the total gene numbers from Autisim
//            Object[] objs = randomizer.nextSample(geneNames, 3900);
//            List<String> sample = InteractionUtilities.convertArrayToList(objs);
//            Map<String, String> geneNameToUniProt = ucscAnalyzer.mapToUniProt(sample);
//            //System.out.println(sample);
//            List<TopicInfo> infos = annotator.topicAnnotate(geneNameToUniProt.values(), 
//                                                            idToTopics, 
//                                                            topicToRatio);
//            for (TopicInfo info : infos) {
//                // Need to check threshold: sorted already
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, System.out);
//    }
//
//    private void outputTestRestults(List<String> hitTopics,
//                                    PrintStream ps) {
//        Map<String, Integer> topicToNumber = InteractionUtilities.countTermUsageInList(hitTopics);
//        for (Iterator<String> it = topicToNumber.keySet().iterator(); it.hasNext();) {
//            String topic = it.next();
//            Integer number = topicToNumber.get(topic);
//            ps.println(topic + "\t" + 
//                               number + "\t" + 
//                               topicToIdNumber.get(topic) + "\t" + 
//                               topicToRatio.get(topic));
//        }
//        ps.println("Total Topics: " + topicToNumber.size());
//    }
//    
//    @Test
//    public void randomTest() throws IOException {
//        List<String> hitTopics = new ArrayList<String>();
//        // Force to use a random seed
//        RandomData randomData = new RandomDataImpl(new JDKRandomGenerator());
//        double threshold = 0.0001d;
//        for (int i = 0; i < 10000; i++) {
//            Set<String> set = new HashSet<String>(idToTopics.keySet());
//            // 909 out of 3767 of total Scott Powers' samples
//            Object[] sampledIds = randomData.nextSample(set, 2000);
//            List<String> sample = new ArrayList<String>();
//            for (Object obj : sampledIds) {
//                sample.add(obj.toString());
//            }
//            //System.out.println(sample);
//            List<TopicInfo> infos = annotator.topicAnnotate(sample, 
//                                                  idToTopics, 
//                                                  topicToRatio);
//            for (TopicInfo info : infos) {
//                // Need to check threshold: sorted already
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, System.out);
//    }
//    
//    /**
//     * A non-redundant random test based on chromosome segments from the provided list of
//     * Segment objects.
//     * @param sampleSegments
//     * @param threshold
//     * @param permutationNumber
//     * @param ps
//     * @throws IOException
//     */
//    public void randomTestBasedOnSegments(List<Segment> sampleSegments,
//                                          double threshold,
//                                          int permutationNumber,
//                                          PrintStream ps) throws IOException {
//        // Want to have chromosome sizes
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        Map<String, Integer> chromosomeSizes = new UCSCDataAnalyzer().loadChromosomeSizes();
//        List<RefSeqInfo> refSeqInfos = ucscAnalyzer.loadHumanRefSeqInfos();
//        Map<String, List<RefSeqInfo>> chrToRefSeqs = ucscAnalyzer.convertRefSeqListToMap(refSeqInfos);
//        List<String> hitTopics = new ArrayList<String>();
//        Random random = new Random();
//        for (int i = 0; i < permutationNumber; i++) {
//            List<Segment> randomSample = generateRandomSegments(sampleSegments,
//                                                                chromosomeSizes,
//                                                                random);
//            Set<String> geneNames = new HashSet<String>();
//            for (Segment segment : randomSample) {
//                List<RefSeqInfo> refSeqList = ucscAnalyzer.grepRefSeqInfo(segment, 
//                                                                          chrToRefSeqs, 
//                                                                          true);
//                //System.out.println(segment + ": " + refSeqList.size());
//                for (RefSeqInfo refSeq : refSeqList)
//                    geneNames.add(refSeq.getGeneName());
//            }
//            Map<String, String> geneNameToUniProt = ucscAnalyzer.mapToUniProt(geneNames);
//            List<TopicInfo> infos = annotator.topicAnnotate(geneNameToUniProt.values(), 
//                                                            idToTopics, 
//                                                            topicToRatio);
//            for (TopicInfo info : infos) {
//                // Need to check threshold: sorted already
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, ps);
//    }
//    
//    /**
//     * A method to do random test based on gene names.
//     * @param geneNames
//     * @param sampleNumber
//     * @param threshold
//     * @param permutation
//     * @param ps
//     * @throws IOException
//     */
//    public void randomTestForGeneNames(Collection<String> geneNames,
//                                       int sampleNumber,
//                                       double threshold,
//                                       int permutation,
//                                       PrintStream ps) throws IOException {
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        RandomData randomizer = new RandomDataImpl(new JDKRandomGenerator());
//        List<String> hitTopics = new ArrayList<String>();
//        for (int i = 0; i < permutation; i++) {
//            Object[] objs = randomizer.nextSample(geneNames, sampleNumber);
//            List<String> sample = InteractionUtilities.convertArrayToList(objs);
//            Map<String, String> geneNameToUniProt = ucscAnalyzer.mapToUniProt(sample);
//            List<TopicInfo> infos = annotator.topicAnnotate(geneNameToUniProt.values(),
//                                                            idToTopics, 
//                                                            topicToRatio);
//            for (TopicInfo info : infos) {
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, ps);
//    }
//    
//    public void redundantRandomTestBasedOnSegments(List<Segment> sampleSegments,
//                                                   Map<String, Integer> segmentFrequence,
//                                                   double threshold,
//                                                   int permutationNumber,
//                                                   PrintStream ps) throws IOException {
//        // Want to have chromosome sizes
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        Map<String, Integer> chromosomeSizes = new UCSCDataAnalyzer().loadChromosomeSizes();
//        List<RefSeqInfo> refSeqInfos = ucscAnalyzer.loadHumanRefSeqInfos();
//        Map<String, List<RefSeqInfo>> chrToRefSeqs = ucscAnalyzer.convertRefSeqListToMap(refSeqInfos);
//        List<String> hitTopics = new ArrayList<String>();
//        Random random = new Random();
//        for (int i = 0; i < permutationNumber; i++) {
//            List<Segment> randomSample = generateRandomSegments(sampleSegments,
//                                                                chromosomeSizes,
//                                                                random);
//            List<String> geneNames = new ArrayList<String>();
//            for (Segment segment : randomSample) {
//                List<RefSeqInfo> refSeqList = ucscAnalyzer.grepRefSeqInfo(segment, 
//                                                                          chrToRefSeqs, 
//                                                                          true);
//                //System.out.println(segment + ": " + refSeqList.size());
//                // depened on how many this segment occurs
//                Integer frequence = segmentFrequence.get(segment.getName());
//                for (int j = 0; j < frequence; j++) {
//                    for (RefSeqInfo refSeq : refSeqList)
//                        geneNames.add(refSeq.getGeneName());
//                }
//            }
//            Map<String, String> geneNameToUniProt = ucscAnalyzer.mapToUniProt(new HashSet<String>(geneNames));
//            List<String> uniProtIds = new ArrayList<String>();
//            for (String geneName : geneNames) {
//                String uniProt = geneNameToUniProt.get(geneName);
//                if (uniProt == null)
//                    continue;
//                uniProtIds.add(uniProt);
//            }
//            List<TopicInfo> infos = annotator.topicAnnotate(uniProtIds, 
//                                                            idToTopics, 
//                                                            topicToRatio);
//            for (TopicInfo info : infos) {
//                // Need to check threshold: sorted already
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, ps);
//    }
//    
//    /**
//     * This method is used to generate a random sample of Segments that have the same
//     * size as in the original samples.
//     */
//    public List<Segment> generateRandomSegments(Collection<Segment> samples,
//                                                Map<String, Integer> chromosomeSizes,
//                                                Random random) {
//        List<Segment> randomSample = new ArrayList<Segment>();
//        for (Segment sampleSegment : samples) {
//            String chromosome = sampleSegment.getChromosome();
//            int totalLength = chromosomeSizes.get(chromosome);
//            int length = sampleSegment.getLength();
//            Segment newSegment = generateRandomSegement(randomSample, 
//                                                        chromosome, 
//                                                        length, 
//                                                        totalLength,
//                                                        random);
//            // Want to keep the same index to be used later
//            newSegment.setName(sampleSegment.getName());
//            randomSample.add(newSegment);
//        }
//        return randomSample;
//    }
//    
//    public List<Segment> generateRandomSegments(List<Segment> samples,
//                                                List<List<String>> geneLists,
//                                                Map<String, Integer> chromosomeSizes,
//                                                UCSCDataAnalyzer ucscAnalyzer,
//                                                Map<String, List<RefSeqInfo>> chrToRefSeqs,
//                                                Random random) {
//        List<Segment> randomSample = new ArrayList<Segment>();
//        int index = 0;
//        List<String> chrList = new ArrayList<String>(chrToRefSeqs.keySet());
//        for (Segment sampleSegment : samples) {
//            String chromosome = sampleSegment.getChromosome();
//            // Use a random chromosom
////            int chrIndex = (int) (chrList.size() * Math.random());
////            String chromosome = chrList.get(chrIndex);
////            System.out.println("Chromosome: " + chromosome);
////            if (!chromosomeSizes.containsKey(chromosome))
////                continue;
//            int totalLength = chromosomeSizes.get(chromosome);
//            int length = sampleSegment.getLength();
//            Segment newSegment = generateRandomSegement(randomSample, 
//                                                        chromosome, 
//                                                        length, 
//                                                        totalLength,
//                                                        geneLists.get(index),
//                                                        ucscAnalyzer,
//                                                        chrToRefSeqs,
//                                                        random);
//            // Want to keep the same index to be used later
//            newSegment.setName(sampleSegment.getName());
//            randomSample.add(newSegment);
//            index ++;
//        }
//        return randomSample;
//    }
//    
//    public List<List<String>> generateRandomGeneSetBasedOnSegmentGenes(List<Segment> samples,
//                                                                       List<List<String>> geneLists,
//                                                                       Map<String, List<String>> chrToGenes) {
//        int index = 0;
//        List<List<String>> rtn = new ArrayList<List<String>>();
//        List<String> chrList = new ArrayList<String>(chrToGenes.keySet());
//        for (Segment sampleSegment : samples) {
//            String chromosome = sampleSegment.getChromosome();
//            //int chrIndex = (int) (Math.random() * chrList.size());
//            //String chromosome = chrList.get(chrIndex);
//            List<String> samplingGenes = chrToGenes.get(chromosome);
//            if (samplingGenes.size() < 100)
//                continue; 
//            List<String> targetGenes = geneLists.get(index);
//            // Get the same number of adjacent genes
//            int geneIndex = (int) (Math.random() * samplingGenes.size());
//            geneIndex -= (targetGenes.size() - 1); // So that we will not get overflow problem
//            if (geneIndex < 0)
//                geneIndex = 0;
//            List<String> list = new ArrayList<String>();
//            for (int i = geneIndex; i < targetGenes.size() + geneIndex; i++)
//                list.add(samplingGenes.get(i));
//            rtn.add(list);
//            index ++;
//        }
//        return rtn;
//    }
//    
//    /**
//     * Use this method to generate a random Segment that will not overlap with any 
//     * other random generated Segments. In the provided Scott Powers' sample, there are
//     * four overlaps among all segments. However, if overlapping is not avoiding in 
//     * this method, 50 or so overlapps will be generated. So for our random testing,
//     * any overlapping is avoided.
//     * @param segments
//     * @param chromosome
//     * @param length
//     * @param totalLength
//     * @return
//     */
//    private Segment generateRandomSegement(List<Segment> segments,
//                                           String chromosome,
//                                           int length,
//                                           int totalLength,
//                                           Random random) {
//        Segment newSegment = new Segment();
//        newSegment.setChromosome(chromosome);
//        // Use this loop to avoid segment overlap
//        boolean isFound = false;
//        // Use one Random only in this method
//        //Random random = new Random();;
//        // To avoid overflow
//        int checkingLength = totalLength - length;
//        //System.out.println("Totallength: " + totalLength + " in " + chromosome);
//        while (!isFound) {
//            isFound = true;
//            int start = random.nextInt(checkingLength);
//            int end = start + length - 1;
//            // Check if it overlaps with any existing segments
//            for (Segment tmp : segments) {
//                if (tmp.overlap(start, end, chromosome)) {
//                    isFound = false;
//                    break;
//                }
//            }
//            if (isFound) {
//                newSegment = new Segment(start,
//                                         end,
//                                         chromosome);
//            }
//        }
//        return newSegment;
//    }
//    
//    private Segment generateRandomSegement(List<Segment> segments,
//                                           String chromosome,
//                                           int length,
//                                           int totalLength,
//                                           List<String> genes,
//                                           UCSCDataAnalyzer ucscAnalyzer,
//                                           Map<String, List<RefSeqInfo>> chrToRefSeqs,
//                                           Random random) {
//        Segment newSegment = new Segment();
//        newSegment.setChromosome(chromosome);
//        // Use this loop to avoid segment overlap
//        boolean isFound = false;
//        // Use one Random only in this method
//        //Random random = new Random();;
//        // To avoid overflow
//        int checkingLength = totalLength - length;
//        Set<String> randomGenes = new HashSet<String>();
//        //System.out.println("Totallength: " + totalLength + " in " + chromosome);
//        while (!isFound) {
//            isFound = true;
//            randomGenes.clear();
//            int start = random.nextInt(checkingLength);
//            int end = start + length - 1;
//            // Check if it overlaps with any existing segments
//            for (Segment tmp : segments) {
//                if (tmp.overlap(start, end, chromosome)) {
//                    isFound = false;
//                    break;
//                }
//            }
//            if (!isFound)
//                continue;
//            newSegment = new Segment(start,
//                                     end,
//                                     chromosome);
//            // Make sure it has the same number of genes
//            List<RefSeqInfo> refSeqList = ucscAnalyzer.grepRefSeqInfo(newSegment, 
//                                                                      chrToRefSeqs, 
//                                                                      false);
//            //System.out.println(segment + ": " + refSeqList.size());
//            for (RefSeqInfo refSeq : refSeqList)
//                randomGenes.add(refSeq.getGeneName());
//            if (randomGenes.size() != genes.size()) {
//                isFound = false;
//            }
//            if (isFound) {
//                break;
//            }
//        }
//        return newSegment;
//    }
//    
//    public void randomTest(Collection<String> totalIds,
//                           int sampleNumber,
//                           double threshold,
//                           int permutationNumber,
//                           PrintStream ps) throws IOException {
//        List<String> hitTopics = new ArrayList<String>();
//        RandomData randomData = new RandomDataImpl(new JDKRandomGenerator());
//        for (int i = 0; i < permutationNumber; i++) {
//            // 909 out from 3767 of total Scott Powers' samples
//            Object[] sampledIds = randomData.nextSample(totalIds, sampleNumber);
//            List<String> sample = new ArrayList<String>();
//            for (Object obj : sampledIds) {
//                sample.add(obj.toString());
//            }
//            //System.out.println(sample);
//            List<TopicInfo> infos = annotator.topicAnnotate(sample, 
//                                                            idToTopics, 
//                                                            topicToRatio);
//            for (TopicInfo info : infos) {
//                // Need to check threshold: sorted already
//                if (info.getPValue() > threshold)
//                    break;
//                hitTopics.add(info.getTopic());
//            }
//        }
//        outputTestRestults(hitTopics, ps);
//    }
//    
//    @Test
//    public void grepPathwayRatiosFromRandomFile() throws IOException {
//        String fileName = R3Constants.RESULT_DIR + "ScottPowers/RandomTestForAllNonRedundant.txt";
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        List<String> titleLines = new ArrayList<String>();
//        String currentTitle = null;
//        int permutationNumber = 0;
//        Map<String, List<String>> pathwayToRatio = new HashMap<String, List<String>>();
//        // This RegExp is used to get the permutation number
//        String regExp = "(\\d+)$";
//        Pattern pattern = Pattern.compile(regExp);
//        while ((line = fu.readLine()) != null) {
//            if (line.equals("------"))
//                break;
//            if (line.length() == 0 ||
//                line.startsWith("Total Topics"))
//                continue;
//            if (line.startsWith("p-value:")) {
//                titleLines.add(line);
//                currentTitle = line;
//                // Need to get permutation number
//                Matcher matcher = pattern.matcher(line);
//                // Have to call find() method first
//                matcher.find();
//                String permutation = matcher.group(1);
//                permutationNumber = Integer.parseInt(permutation);
//            }
//            else {
//                String[] tokens = line.split("\t");
//                // Translation(R)   7   119 0.018772677078403535
//                double ratio = Double.parseDouble(tokens[1]) / permutationNumber;
//                List<String> list = pathwayToRatio.get(tokens[0]);
//                if (list == null) {
//                    list = new ArrayList<String>();
//                    pathwayToRatio.put(tokens[0], list);
//                }
//                list.add(currentTitle + ", " + ratio);
//            }
//        }
//        fu.close();
//        System.out.print("Pathway\t");
//        for (Iterator<String> it = titleLines.iterator(); it.hasNext();) {
//            line = it.next();
//            System.out.print(line);
//            if (it.hasNext())
//                System.out.print("\t");
//        }
//        // Close this line
//        System.out.println();
//        for (Iterator<String> it = pathwayToRatio.keySet().iterator(); it.hasNext();) {
//            String pathway = it.next();
//            List<String> values = pathwayToRatio.get(pathway);
//            System.out.print(pathway + "\t");
//            for (int i = 0; i < titleLines.size(); i++) {
//                String title = titleLines.get(i);
//                // Check values
//                for (String value : values) {
//                    if (value.startsWith(title)) {
//                        // Get the number
//                        int index = value.lastIndexOf(',');
//                        System.out.print(value.substring(index + 1));
//                        break;
//                    }
//                }
//                System.out.print("\t");
//            }
//            System.out.println();
//        }
//    }
//    
//    @Test
//    public void processRandomFile() throws IOException {
//        //String fileName = "/Users/wgm/Documents/caBIG_R3/xiaoyue/RandomTestForDenovoGenes.txt";
//        //String fileName = R3Constants.RESULT_DIR  + "RandomPathwayHitTest_10000.txt";
//        String fileName = R3Constants.RESULT_DIR + "ScottPowers/RandomTestForAllRedundantBasedOnSegments.txt";
//        int permutationNumber = 1000;
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        int total = 0;
//        String pValueLine = null;
//        while ((line = fu.readLine()) != null) {
//            if (line.equals("------"))
//                break;
//            if (line.length() == 0 ||
//                line.startsWith("Total Topics"))
//                continue;
//            if (line.startsWith("p-value:")) {
//                if (pValueLine != null) {
//                    System.out.println((pValueLine + ": " + total + ": " + total / (double)permutationNumber).replaceAll(":", ""));
//                    total = 0;
//                }
//                pValueLine = line;
//            }
//            else {
//                String[] tokens = line.split("\t");
//                total += Integer.parseInt(tokens[1]);
//            }
//        }
//        fu.close();
//        // Need to print out the last line
//        System.out.println((pValueLine + ": " + total + ": " + total / (double)permutationNumber).replaceAll(":", ""));
//    }
}
