/*
 * Created on Jun 6, 2013
 *
 */
package org.reactome.cancer;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.imageio.ImageIO;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.UCSCDataAnalyzer;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.fi.SurvivalAnalysisResult;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to analyze CNV based breast cancer athracycline resistence.
 * @author gwu
 *
 */
public class BreastCancerAnthracyclineAnalyzer {
    
    private final String DIR_NAME = "datasets/AntResBreastCancer/";
    private final String TCGA_DIR_NAME = "datasets/TCGA/BRCA/06231013/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2013052300.0.0/";
    private final FileUtility fu = new FileUtility();
    
    public BreastCancerAnthracyclineAnalyzer() {
    }
    
    /**
     * Generate a map from chosen genes to types of CNVs (deletation or amplification)
     * @throws Exception
     */
    @Test
    public void generateGeneToCNVTypeMap() throws Exception {
        Map<String, Set<String>> cnvToGenes = loadCNVToGenes();
        Map<String, Set<String>> geneToCNVs = InteractionUtilities.switchKeyValues(cnvToGenes);
        Set<String> chosenGenes = fu.loadInteractions(DIR_NAME + "ChosenGeneList_101613.txt");
        for (String gene : chosenGenes) {
            Set<String> cnvs = geneToCNVs.get(gene);
            Set<String> types = new HashSet<String>();
            for (String cnv : cnvs)
                types.add(cnv.substring(0, 1));
            if (types.size() > 1) {
//                System.err.println("Gene has more than one type: " + gene + " (" + cnvs + ")");
                // The following fours will be detected: A159, A13, D69, A50
                // Choose A since it is used in GroupB
                System.out.println(gene + "\tA");
                continue;
            }
            System.out.println(gene + "\t" + types.iterator().next());
        }
    }
    
    @Test
    public void generateGeneList() throws IOException {
        String fileName = DIR_NAME + "ChosenGenes_101613.txt";
        Set<String> genes = new HashSet<String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        fu.close();
        String[] tokens = line.split(", ");
        for (String token : tokens)
            genes.add(token);
        System.out.println("Total genes: " + genes.size());
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        fileName = DIR_NAME + "ChosenGeneList_101613.txt";
        fu.saveCollection(geneList, fileName);
    }
    
    @Test
    public void chooseTargetGenesForFutherStudies() throws Exception {
        Set<String> selectedGenes = new HashSet<String>();
        
        // Load Group A genes
        System.out.println("Choosing genes from Group A:");
        String fileName = DIR_NAME + "GenesInA61_A86_A132_A162_A203_092713.txt";
        Set<String> groupAGenes = fu.loadInteractions(fileName);
        System.out.println("Total genes in group A: " + groupAGenes.size());
        // Group A genes in the FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> groupAGenesInFI = InteractionUtilities.getShared(fiGenes, groupAGenes);
        System.out.println("Genes in FI: " + groupAGenesInFI.size());
        // Genes that interact each other
        Set<String> fisInGroupA = InteractionUtilities.getFIs(groupAGenes, fis);
        Set<String> groupAFIGenes = InteractionUtilities.grepIDsFromInteractions(fisInGroupA);
        System.out.println("Interacting genes: " + groupAFIGenes.size());
        // Do a pathway analysis in order to pick up pathway related genes
        PathwayBasedAnnotator pathwayAnnotator = new PathwayBasedAnnotator();
        pathwayAnnotator.setFDRThreshold(0.1d); // This is based on a pre-checked results so that four pathways can be picked.
        List<GeneSetAnnotation> annotations = pathwayAnnotator.annotateGenesWithFDR(groupAGenes,
                                                                                    AnnotationType.Pathway);
        CancerAnalysisUtilitites.printOutAnnotations(annotations);
        // The first two pathways are related to Olfatory receptor gene cluster.
        List<String> genes1 = annotations.get(0).getHitIds();
        List<String> genes2 = annotations.get(1).getHitIds();
        Set<String> shared = InteractionUtilities.getShared(genes1, genes2);
        genes1.removeAll(shared);
        genes2.removeAll(shared);
        Set<String> groupAGenesFromPathways = new HashSet<String>();
        groupAGenesFromPathways.addAll(genes1);
        groupAGenesFromPathways.addAll(genes2);
        groupAGenesFromPathways.addAll(annotations.get(2).getHitIds());
        groupAGenesFromPathways.addAll(annotations.get(3).getHitIds());
        System.out.println("Total pathway genes without ORs: " + groupAGenesFromPathways.size());
        // Add pathway genes and FI genes together
        Set<String> selectedGroupAGenes = new HashSet<String>();
        selectedGroupAGenes.addAll(groupAGenesFromPathways);
        selectedGroupAGenes.addAll(groupAFIGenes);
        // Get half of ORs genes
        Set<String> halfShared = MathUtilities.randomSampling(shared, shared.size() / 2);
        selectedGroupAGenes.addAll(halfShared);
        // We should get 131 genes from Group A
        System.out.println("Total selected Group A genes: " + selectedGroupAGenes.size());
        selectedGenes.addAll(selectedGroupAGenes);
        List<String> geneList = new ArrayList<String>(selectedGroupAGenes);
        Collections.sort(geneList);
        System.out.println(geneList);
        
        // Get all genes from Group B
        System.out.println("\nChoosing genes from Group B:");
        fileName = DIR_NAME + "GenesInA13_D18_D47_092713.txt";
        Set<String> groupBGenes = fu.loadInteractions(fileName);
        System.out.println("Group B genes: " + groupBGenes.size());
        Set<String> groupBGenesInFI = InteractionUtilities.getShared(fiGenes, groupBGenes);
        System.out.println("Genes in FIs: " + groupBGenesInFI.size());
        Set<String> selectedGroupBGenes = new HashSet<String>();
        // Choose all genes in the FI network: 30 genes in total
        selectedGroupBGenes.addAll(groupBGenesInFI);
        // Pick up some random genes from each individual fragments
        Map<String, Set<String>> cnvToGenes = loadCNVToGenes();
        String[] cnvs = new String[] {"A13", "D18", "D47"};
        for (String cnv : cnvs) {
            Set<String> genes = cnvToGenes.get(cnv);
            System.out.println(cnv + "\t" + genes.size());
            // Get 1/2 genes from each CNV
            selectedGroupBGenes.addAll(MathUtilities.randomSampling(genes, genes.size() / 2));
        }
        System.out.println("Total selected Group B genes: " + selectedGroupBGenes.size());
        geneList = new ArrayList<String>(selectedGroupBGenes);
        Collections.sort(geneList);
        System.out.println(geneList);
        selectedGenes.addAll(selectedGroupBGenes);
        
        // Group C doesn't have a lot of CNV+ samples. 69 genes in total in this group. Only choose 1/4 genes.
        System.out.println("\nChoosing genes from Group C:");
        fileName = DIR_NAME + "GenesInA59_A73_D29_092713.txt";
        Set<String> groupCGenes = fu.loadInteractions(fileName);
        System.out.println("Group C genes: " + groupCGenes.size());
        cnvs = new String[] {"A59", "A73", "D29"};
        Set<String> selectedGroupCGenes = new HashSet<String>();
        for (String cnv : cnvs) {
            Set<String> genes = cnvToGenes.get(cnv);
            System.out.println(cnv + "\t" + genes.size());
            // Get 1/4 genes from each CNV
            selectedGroupCGenes.addAll(MathUtilities.randomSampling(genes, genes.size() / 4));
        }
        System.out.println("Total selected Group C genes: " + selectedGroupCGenes.size());
        geneList = new ArrayList<String>(selectedGroupCGenes);
        Collections.sort(geneList);
        System.out.println(geneList);
        selectedGenes.addAll(selectedGroupCGenes);
        
        // Group D doesn't have enough CNV+ samples. Will not consider
        System.out.println("\nChoosing genes from Group D: No genes are chosen because Group D doesn't have enough CNV+ samples!");
        
        // Group E: Samples are not enough. It has 246 genes in total. Choose 1/4 genes.
        System.out.println("\nChoosing genes from Group E:");
        fileName = DIR_NAME + "GenesInD1_D30_D45_092713.txt";
        Set<String> groupEGenes = fu.loadInteractions(fileName);
        System.out.println("Total Group E genes: " + groupEGenes.size());
        Set<String> selectedGroupEGenes = new HashSet<String>();
        // Get all genes in enriched pathways
        pathwayAnnotator.setFDRThreshold(0.25d);
        annotations = pathwayAnnotator.annotateGenesWithFDR(groupEGenes,
                                                            AnnotationType.Pathway);
        for (GeneSetAnnotation annotation : annotations) {
            System.out.println(annotation.getTopic() + "\t" + annotation.getHitNumber());
            selectedGroupEGenes.addAll(annotation.getHitIds());
        }
//        cnvs = new String[]{"D1", 
////                            "D30", 
//                            "D45"}; // D1 is contained by D30. Choose genes from D1 only
//        for (String cnv : cnvs) {
//            Set<String> genes = cnvToGenes.get(cnv);
//            System.out.println(cnv + "\t" + genes.size());
//            // Get 1/3 genes from each CNV
//            selectedGroupEGenes.addAll(MathUtilities.randomSampling(genes, genes.size() / 4));
//        }
        // Only 1 gene in D45. Select it
        Set<String> d45Genes = cnvToGenes.get("D45");
        System.out.println("D45 genes: " + d45Genes);
        selectedGroupEGenes.addAll(cnvToGenes.get("D45")); 
        System.out.println("Total selected Group E genes: " + selectedGroupEGenes.size());
        geneList = new ArrayList<String>(selectedGroupEGenes);
        Collections.sort(geneList);
        System.out.println(geneList);
        selectedGenes.addAll(selectedGroupEGenes);
        
        // Group F has 21 genes only. Select all of them.
        System.out.println("\nChoosing genes from Group F:");
        fileName = DIR_NAME + "GenesInA89_A110_092713.txt";
        Set<String> groupFGenes = fu.loadInteractions(fileName);
        System.out.println("Group F genes: " + groupFGenes.size());
        Set<String> selectedGroupFGenes = new HashSet<String>(groupFGenes);
        System.out.println("Total selected Gropu F genes: " + selectedGroupFGenes.size());
        geneList = new ArrayList<String>(selectedGroupFGenes);
        Collections.sort(geneList);
        System.out.println(geneList);
        selectedGenes.addAll(selectedGroupFGenes);
        
        System.out.println("\n Total selected genes: " + selectedGenes.size());
        geneList = new ArrayList<String>(selectedGenes);
        Collections.sort(geneList);
        for (String gene : geneList)
            System.out.println(gene);
    }
    
    @Test
    public void checkOlfactoryAndGPCRGenesForGroupA() {
        String olfGenes = "OR2T2, OR2T1, OR2T4, OR2T3, OR2T5, OR2T6, OR2T8, OR2W3, OR2AK2, OR2L8, OR2L3, OR2L13, OR2L2, OR2T33, OR2T34, OR2M7, OR2M2, OR2M3, OR2M4, OR2M5, OR2T29, OR2T27, OR2T10, OR2T11, OR2T12, OR2B11, PDC, OR13G1, OR2C3, OR1C1, OR6F1, OR2G3, OR2G2, OR2G6, OR14C36, OR14A16, OR14I1, OR11L1";
        String gpcrGenes = "GNG4, OR2T2, OR2T1, OR2T4, OR2T3, OR2T5, OR2T6, OR2T8, RGS1, RGS2, AVPR1B, RGS7, RGS8, KISS1, OR2W5, OR2W3, OR2AK2, OR2L8, OR2L5, OR2L3, OR2L13, OR2L2, OR2T33, OR2T35, OR2T34, OR2M7, OR2M2, OR2M3, OR2M4, OR2M5, OR2T29, OR2T27, XCL1, XCL2, OR2T10, OR2T11, OR2T12, CHRM3, OR2B11, ADORA1, OR2C3, OR1C1, OR6F1, OR2G3, OR2G2, OR2G6, AGT, OPN3, OR14C36, OBSCN, OR14A16, AKT3, OR14I1, RGS18, RGS21";
        String[] tokens = olfGenes.split(", ");
        List<String> olfGeneList = Arrays.asList(tokens);
        tokens = gpcrGenes.split(", ");
        List<String> gpcrGeneList = Arrays.asList(tokens);
        System.out.println("Olfactory genes: " + olfGeneList.size());
        System.out.println("GPCR genes: " + gpcrGeneList.size());
        Set<String> shared = InteractionUtilities.getShared(gpcrGeneList, olfGeneList);
        System.out.println("Shared: " + shared.size());
        List<String> sharedList = new ArrayList<String>(shared);
        Collections.sort(sharedList);
        System.out.println("Shared genes: " + StringUtils.join(sharedList.iterator(), ", "));
    }
    
    @Test
    public void generateCNVGroupsToGenes() throws Exception {
        Map<String, Set<String>> cnvToGenes = loadCNVToGenes();
        // The following are pre-selected CNVs groups
        String[] cnvGroups = new String[] {
                "A61, A86, A132, A162, A203",
                "A13, D18, D47",
                "A59, A73, D29",
                "D43, D52, D53",
                "D1, D30, D45",
                "A89, A110" // A198 genes have not checked
        };
        PathwayBasedAnnotator pathwayAnnotator = new PathwayBasedAnnotator();
        pathwayAnnotator.setFDRThreshold(0.5d);
        for (String cnvGroup : cnvGroups) {
            String[] cnvs = cnvGroup.split(", ");
            Set<String> genes = new HashSet<String>();
            for (String cnv : cnvs) {
                Set<String> cnvGenes = cnvToGenes.get(cnv);
                if (cnvGenes == null || cnvGenes.size() == 0) {
                    System.out.println(cnv + " has no genes!");
                    continue;
                }
                genes.addAll(cnvToGenes.get(cnv));
            }
            String title = StringUtils.join(cnvs, "_");
            String fileName = DIR_NAME + "GenesIn" + title + "_092713.txt";
            fu.saveCollection(genes, fileName);
            // Want to do a pathway enrichment analysis
            System.out.println("Pathawys for " + title);
            List<GeneSetAnnotation> annotations = pathwayAnnotator.annotateGenesWithFDR(genes,
                                                                                       AnnotationType.Pathway);
            CancerAnalysisUtilitites.printOutAnnotations(annotations);
            System.out.println();
        }
    }
    
    @Test
    public void drawCNVFragments() throws IOException {
        String text = "A13, A59, A61, A73, A86, A132, A162, A203, D18, D29, D43, D45, D47, D52, D53, D1, D30, A89, A198, A110";
        String[] tokens = text.split(", ");
        List<String> cnvs = Arrays.asList(tokens);
        System.out.println("Total selected CNVs: " + cnvs.size());
        Map<String, Integer> chromToLength = loadChromToLength();
        Map<String, String> cnvToLocation = loadFragmentToLocations();
        for (String chrom : chromToLength.keySet()) {
            Integer length = chromToLength.get(chrom);
            System.out.println(chrom + ": " + length);
            List<String> cnvsToBeDrawn = new ArrayList<String>();
            for (String cnv : cnvs) {
                String location = cnvToLocation.get(cnv);
                int index = location.indexOf(":");
                String chrom1 = location.substring(0, index);
                if (chrom1.equals(chrom)) {
                    System.out.println(cnv + ": " + location);
                    cnvsToBeDrawn.add(cnv);
                }
            }
            System.out.println();
            if (cnvsToBeDrawn.size() == 0)
                continue;
            int width = 1200;
            int height = 600;
            int buffer = 20;
            BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
            Graphics2D g2 = (Graphics2D) image.getGraphics();
            g2.fillRect(0, 0, width, height);
            drawChromosome(chrom, 
                           length,
                           width,
                           height,
                           buffer,
                           g2);
            double legnthToCoor = (width - 2 * buffer) / (double) length;
            int step = 20;
            int y1 = height / 2 - 40; // For amplification
            int y2 = height / 2 + 40; // For deletion 
            FontMetrics fm = g2.getFontMetrics();
            for (String cnv : cnvsToBeDrawn) {
                String location = cnvToLocation.get(cnv);
                int index = location.indexOf(":");
                location = location.substring(index + 1);
                index = location.indexOf("-");
                int loc1 = new Integer(location.substring(0, index));
                int loc2 = new Integer(location.substring(index + 1));
                int x1 = (int) (loc1 * legnthToCoor) + buffer;
                int x2 = (int) (loc2 * legnthToCoor) + buffer;
                if (cnv.startsWith("A")) {
                    g2.drawLine(x1, y1, x2, y1);
                    // Draw labels
                    int txtWidth = fm.stringWidth(cnv);
                    int txtX = (x1 + x2 - txtWidth) / 2;
                    g2.drawString(cnv, txtX, y1 - 2);
                    y1 -= step;
                }
                else {
                    g2.drawLine(x1, y2, x2, y2);
                    // Draw labels
                    int txtWidth = fm.stringWidth(cnv);
                    int txtX = (x1 + x2 - txtWidth) / 2;
                    g2.drawString(cnv, txtX, y2 + 14);
                    y2 += step;
                }
            }
            ImageIO.write(image,
                          "png",
                          new File(DIR_NAME + "SelectedCNVs_" + chrom + "_092713.png"));
        }
    }
    
    private Map<String, String> loadFragmentToLocations() throws IOException {
        Map<String, String> fragmentToLocation = new HashMap<String, String>();
        String scoreFile = DIR_NAME + "NTcoresCoxPh_092313.txt";
        fu.setInput(scoreFile);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // These two cores should be removed because they are duplicates of A13
            if (tokens[0].equals("A50") ||
                tokens[0].equals("A159"))
                continue;
            fragmentToLocation.put(tokens[0],
                                   tokens[2]);
        }
        fu.close();
        return fragmentToLocation;
    }
    
    private void drawChromosome(String chr,
                                Integer length,
                                int width,
                                int height,
                                int buffer,
                                Graphics2D g2) {
        // Draw the central line
        int x1 = buffer;
        int y1 = height / 2;
        int x2 = width - buffer;
        int y2 = y1;
        g2.setPaint(Color.BLACK);
        g2.drawLine(x1, y1, x2, y2);
        // Draw 10 marks
        int marks = 10;
        int step = (width - 2 * buffer) / marks;
        int lengthStep = length / marks;
        for (int i = 0; i < marks + 1; i++) { // +1 for the last tick
            g2.drawLine(x1, y2 + 5, x1, y2 - 5); // Double ticks
            x1 += step;
            // Draw length
            if (i < marks - 1) {
                String tick = lengthStep * (i + 1) + "";
                int tickLength = g2.getFontMetrics().stringWidth(tick);
                x2 = x1 - tickLength / 2;
                g2.drawString(tick, 
                              x2, 
                              y2 + 20);
            }
        }
        // Draw a title
        int index = chr.indexOf("r"); // for 1 or 11.
        String title = "Chromosome " + chr.substring(index + 1);
        int titleWidth = g2.getFontMetrics().stringWidth(title);
        x1 = (width - titleWidth) / 2;
        y1 = 25;
        g2.drawString(title, x1, y1);
    }
    
    @Test
    public void drawCNVsInChromosomes() throws IOException {
        String line;
        Map<String, Integer> chrToLength = loadChromToLength();
//        for (String chr : chrToLength.keySet())
//            System.out.println(chr + "\t" + chrToLength.get(chr));
        Map<String, String> fragmentToLocation = new HashMap<String, String>();
        Map<String, Double> fragmentToScore = new HashMap<String, Double>();
        String scoreFile = DIR_NAME + "NTcoresCoxPh_092313.txt";
        fu.setInput(scoreFile);
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // These two cores should be removed because they are duplicates of A13
            if (tokens[0].equals("A50") ||
                tokens[0].equals("A159"))
                continue;
            fragmentToLocation.put(tokens[0],
                                   tokens[2]);
            fragmentToScore.put(tokens[0],
                                new Double(tokens[3]));
        }
        fu.close();
        // Check the duplication
//        System.out.println("Total fragments: " + fragmentToLocation.size());
//        List<String> locationSet = new ArrayList<String>(fragmentToLocation.values());
//        Map<String, Integer> locationToCount = InteractionUtilities.countTermUsageInList(locationSet);
//        Set<String> duplicated = new HashSet<String>();
//        for (String location : locationToCount.keySet()) {
//            Integer count = locationToCount.get(location);
//            if (count > 1)
//                System.out.println("Duplicated: " + location);
//        }
        List<String> chromList = new ArrayList<String>(chrToLength.keySet());
        Collections.sort(chromList);
        int width = 1200;
        int height = 600;
        int buffer = 20;
        for (String chr : chromList) {
            if (!chr.toLowerCase().startsWith("chr8"))
                continue;
            Integer length = chrToLength.get(chr);
            System.out.println("Chr: " + length);
            BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
            Graphics2D g2 = (Graphics2D) image.getGraphics();
            g2.fillRect(0,
                        0, 
                        image.getWidth(), 
                        image.getHeight());
            drawChromosome(chr, length, width, height, buffer, g2);
            // Draw the actual values
            Map<int[], Double> ampLocationToScore = new HashMap<int[], Double>();
            Map<int[], Double> delLocationToScore = new HashMap<int[], Double>();
            for (String fragment : fragmentToLocation.keySet()) {
                String location = fragmentToLocation.get(fragment);
                int index = location.indexOf(":");
                String coreChr = location.substring(0, index);
                if (coreChr.equals(chr)) {
                    System.out.println(fragment + "\t" + 
                                       location.substring(index + 1) + "\t" +
                                       fragmentToScore.get(fragment));
                    int[] coors = new int[2];
                    location = location.substring(index + 1);
                    String[] coorText = location.split("-");
                    coors[0] = new Integer(coorText[0]);
                    coors[1] = new Integer(coorText[1]);
                    if (fragment.startsWith("A"))
                        ampLocationToScore.put(coors, 
                                               fragmentToScore.get(fragment));
                    else
                        delLocationToScore.put(coors, 
                                               fragmentToScore.get(fragment));
                }
            }
            // Map from length to coordinate
            double lengthToCoor = (double) (width - 2 * buffer) / length;
            double ampMax = 0.0d;
            int currentX = buffer;
            double scoreToHeight = 7.5d;
            SummaryStatistics stat = new SummaryStatistics();
            // Use some transparent colors
            g2.setPaint(new Color(0.0f, 0.0f, 0.0f, 0.5f));
            String fileName = DIR_NAME + "NTCoresIn" + chr + ".txt";
            fu.setOutput(fileName);
            fu.printLine("browser position " + chr.toLowerCase() + ":1-" + chrToLength.get(chr));
            fu.printLine("browser hide all");
            fu.printLine("track name=BreastCancerCGH(Amp) description=\"Breast Cancer CGH\" visibility=dense useScore=1");
            // Draw amplication scores
            int start = 1;
            Double prevValue = null;
            int count = 1;
            for (int i = 1; i < length; i++) {
                double totalScore = 0.0d;
                for (int[] location : ampLocationToScore.keySet()) {
                    if (i >= location[0] && i <= location[1]) {
                        totalScore += ampLocationToScore.get(location);
                    }
                }
                if (prevValue == null)
                    prevValue = totalScore;
                if (totalScore != prevValue) {
                    int score = (int) (prevValue * 25);
                    if (score > 0) {
                        fu.printLine(chr.toLowerCase() + "\tCore\tCNV\t" +
                                start + "\t" + 
                                i + "\t" + 
                                score + "\t.\t.\tA_" + count++);
                    }
                    start = i;
                    prevValue = totalScore;
                }
//                x1 = (int)(lengthToCoor * i) + buffer;
//                if (x1 == currentX) {
//                    stat.addValue(totalScore);
//                }
//                else {
//                    y1 = (int)(stat.getMean() * scoreToHeight);
//                    g2.drawLine(currentX, y2 - y1, currentX, y2);
//                    currentX = x1;
//                    stat.clear();
//                }
                if (totalScore > ampMax)
                    ampMax = totalScore;
            }
            System.out.println("Max amplication: " + ampMax);
            // Draw deletion scores
            stat.clear();
            currentX = buffer;
            double delMax = 0.0d;
            fu.printLine("track name=BreastCancerCGH(Del) description=\"Breast Cancer CGH\" visibility=dense useScore=1");
            count = 1;
            prevValue = null;
            start = 1;
            for (int i = 1; i < length; i++) {
                double totalScore = 0.0d;
                for (int[] location : delLocationToScore.keySet()) {
                    if (i >= location[0] && i <= location[1]) {
                        totalScore += delLocationToScore.get(location);
                    }
                }
                if (prevValue == null)
                    prevValue = totalScore;
                if (totalScore != prevValue) {
                    int score = (int) (prevValue * 25);
                    if (score > 0) {
                        fu.printLine(chr.toLowerCase() + "\tCore\tCNV\t" +
                                start + "\t" + 
                                i + "\t" + 
                                score + "\t.\t.\tD_" + count++);
                    }
                    start = i;
                    prevValue = totalScore;
                }
//                x1 = (int)(lengthToCoor * i) + buffer;
//                if (x1 == currentX) {
//                    stat.addValue(totalScore);
//                }
//                else {
//                    y1 = (int)(stat.getMean() * scoreToHeight);
//                    g2.drawLine(currentX, y2 + y1, currentX, y2);
//                    currentX = x1;
//                    stat.clear();
//                }
                if (totalScore > delMax)
                    delMax = totalScore;
            }
            System.out.println("Max deletion: " + delMax);
            // Draw vertical lines
            g2.setPaint(Color.BLACK);
            int y2 = height / 2;
            int ampMaxY = (int)(y2 - scoreToHeight * ampMax);
            int delMaxY = (int)(y2 + scoreToHeight * delMax);
            g2.drawLine(buffer, 
                        ampMaxY, 
                        buffer, 
                        delMaxY);
            // Draw vertical ticks and labels
            g2.drawLine(buffer, 
                        ampMaxY, 
                        buffer + 5, 
                        ampMaxY);
            g2.drawLine(buffer,
                        delMaxY,
                        buffer + 5,
                        delMaxY);
            // Draw two labels
            g2.drawString(String.format("%.2f", ampMax),
                          buffer + 7, 
                          ampMaxY + 5);
            g2.drawString(String.format("%.2f", delMax),
                          buffer + 7,
                          delMaxY + 5);
//            ImageIO.write(image, 
//                          "png", 
//                          new File(DIR_NAME + "NTCoresIn" + chr + ".png"));
            fu.close();
//            break;
        }
    }

    private Map<String, Integer> loadChromToLength() throws IOException {
        String chromoFile = DIR_NAME + "ChromosomeLengths_092313.txt";
        fu.setInput(chromoFile);
        String line = fu.readLine();
        Map<String, Integer> chrToLength = new HashMap<String, Integer>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String length = tokens[1].replaceAll(",", "");
            chrToLength.put("Chr" + tokens[0],
                            new Integer(length));
        }
        fu.close();
        return chrToLength;
    }
    
    @Test
    public void attachCoxPHScoresToNTCores() throws IOException {
        String coxphFile = DIR_NAME + "CoxPHByFourGroups_073113.txt";
        Map<String, Double> coreToScore = new HashMap<String, Double>();
        fu.setInput(coxphFile);
        String line = fu.readLine(); // Escape the first two lines
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            coreToScore.put(tokens[0],
                            -Math.log(new Double(tokens[5]))); // Use e as the base for larger separations of scores.
        }
        fu.close();
//        for (String core : coreToScore.keySet())
//            System.out.println(core + "\t" + coreToScore.get(core));
        String coreFile = DIR_NAME + "NTcores.txt";
        String outFile = DIR_NAME + "NTcoresCoxPh_092313.txt";
        fu.setInput(coreFile);
        fu.setOutput(outFile);
        line = fu.readLine();
        fu.printLine(line + "\tCoxPHScore");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double coxphScore = coreToScore.get(tokens[0]);
            fu.printLine(line + "\t" + coxphScore);
        }
        fu.close();
    }
    
    private Map<String, Set<String>> loadCNVToGenes() throws IOException {
        String fileName = DIR_NAME + "GenesInNTCores_080113.txt";
        Map<String, Set<String>> cnvToGenes = new HashMap<String, Set<String>>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String[] genes = tokens[2].split(", ");
            Set<String> geneSet = new HashSet<String>();
            for (String gene : genes)
                geneSet.add(gene);
            cnvToGenes.put(tokens[0], geneSet);
        }
        fu.close();
        return cnvToGenes;
    }
    
    @Test
    public void geneOverlappingAnalysis() throws Exception {
        Map<String, List<String[]>> geneToCoordinates = new UCSCDataAnalyzer().loadGeneNameToCoordinates();
        int totalGene = geneToCoordinates.size();
        System.out.println("Total genes: " + totalGene);
        Map<String, Set<String>> cnvToGenes = new HashMap<String, Set<String>>();
        String fileName = DIR_NAME + "GenesInNTCores_080113.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String[] genes = tokens[2].split(", ");
            Set<String> geneSet = new HashSet<String>();
            for (String gene : genes)
                geneSet.add(gene);
            cnvToGenes.put(tokens[0], geneSet);
        }
        fu.close();
//        for (String cnv : cnvToGenes.keySet()) {
//            Set<String> genes = cnvToGenes.get(cnv);
//            System.out.println(cnv + "\t" + genes.size());
//        }
        // Want to print out for a list of CNVs tha have been subject to Kaplan-Meier survival analysis
        fileName = DIR_NAME + "KaplanMeierSurvivalForFourGroups_073113.txt";
        fu.setInput(fileName);
        List<String> selectedCNVs = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                continue;
            if (line.matches("(A|D)(\\d+)")) {
                selectedCNVs.add(line);
            }
        }
        fu.close();
        System.out.println("Total selected CNVs: " + selectedCNVs.size() + ": " + selectedCNVs);
        // Do hypergeometric text based overlapping analysis
        fileName = DIR_NAME + "CNVFragmentGenesOverlap_080113.txt";
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        for (String cnv : selectedCNVs)
            builder.append("\t").append(cnv);
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String cnv : selectedCNVs) {
            builder.append(cnv);
            Set<String> genes = cnvToGenes.get(cnv);
            for (String cnv1 : selectedCNVs) {
                Set<String> genes1 = cnvToGenes.get(cnv1);
                if (cnv.equals(cnv1)) {
                    builder.append("\tNA");
                }
                else {
                    Set<String> shared = InteractionUtilities.getShared(genes, genes1);
                    if (shared.size() == 0) {
                        builder.append("\t").append(shared.size());
                    }
                    else {
                        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGene,
                                                                                    genes.size(),
                                                                                    genes1.size(), 
                                                                                    shared.size());
                        builder.append("\t").append(shared.size()).append("(").append(String.format("%.2e", pvalue)).append(")");
                    }
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkGenesInCNVFragments() throws IOException {
        Map<String, List<String[]>> geneToCoordinates = new UCSCDataAnalyzer().loadGeneNameToCoordinates();
        String srcFileName = DIR_NAME + "NTcores.txt";
        fu.setInput(srcFileName);
        String line = fu.readLine();
        System.out.println("CNVCore\tNumberOfGenes\tGenes");
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
//            if (!tokens[0].equals("A203"))
//                continue;
            String[] location = tokens[2].split("(:|-)");
            // Only need to lower the first letter since we need to keep X,Y as upper
            // case
            String chr = location[0];
            chr = chr.substring(0, 1).toLowerCase() + chr.substring(1);
            int start = Integer.parseInt(location[1]);
            int end = Integer.parseInt(location[2]);
            Set<String> genes = new HashSet<String>();
            for (String gene : geneToCoordinates.keySet()) {
                List<String[]> list = geneToCoordinates.get(gene);
                for (String[] coord : list) {
                    if (!coord[0].equals(chr))
                        continue;
                    int start0 = Integer.parseInt(coord[1]);
                    int end0 = Integer.parseInt(coord[2]);
                    // Check if there any overlap
                    if (end > start0 & start < end0) {
                        genes.add(gene);
                    }
                }
            }
            List<String> geneList = new ArrayList<String>(genes);
            Collections.sort(geneList);
            System.out.println(tokens[0] + "\t" + 
                               geneList.size() + "\t" +
                               InteractionUtilities.joinStringElements(", ", geneList));
        }
    }
    
    @Test
    public void doTCGABRCAMCLModuleSurvivalAnalysis() throws Exception {
        String scoreFile = DIR_NAME + "TCGA_BRCA_SampleToMCLModuleToValue_AllSamples_070913.txt";
        String clinFile = DIR_NAME + "TCGA_BRCA.clin.merged.picked.transformed.txt";
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        SurvivalAnalysisResult result = survivalHelper.doSurvivalAnalysis(new File(scoreFile), 
                                                                          new File(clinFile), 
                                                                          "coxph", 
                                                                          null, 
                                                                          "tcga_",
                                                                          false);
        String output = CancerAnalysisUtilitites.attachFDRsToSurvivalResults(result.getOutput());
        System.out.println(output);
    }
    
    /**
     * Generate a test data set using the TCGA breast cancer data set using MCL modules
     * selected from CSHL CGH breast cancer data via superpc.
     * @throws Exception
     */
    @Test
    public void generateSampleToScoresForSelectedModules() throws Exception {
        // For all samples
//        String selectedModules = "40, 6, 46, 14, 45, 47, 41";
        //String fileName = DIR_NAME + "CoxPHMCLModulesLongerDRFS_AllSamples_070313.txt";
        String fileName = DIR_NAME + "CoxPHMCLModulesLongerDRFS_AllSamples_070913.txt";
        Map<String, Integer> geneToModule = selectMCLModules(null, 
                                                             fileName, 
                                                             null, 
                                                             null);
        Set<Integer> moduleIndices = new HashSet<Integer>(geneToModule.values());
        List<Integer> moduleList = new ArrayList<Integer>(moduleIndices);
        Collections.sort(moduleList);
        List<String> moduleNames = new ArrayList<String>(moduleList.size());
        for (Integer module : moduleList)
            moduleNames.add("Module" + module);
        System.out.println("Modules: " + moduleNames);
        // Generate back to modules
        List<Set<String>> modules = new ArrayList<Set<String>>();
        // Inititalize a list of empty modules
        for (int i = 0; i < moduleList.size(); i++)
            modules.add(new HashSet<String>());
        for (String gene : geneToModule.keySet()) {
            Integer moduleIndex = geneToModule.get(gene);
            Integer index = moduleList.indexOf(moduleIndex);
            Set<String> module = modules.get(index);
            module.add(gene);
        }
        // Load TCGA Breast cancer data set
        String tcgaFileName = TCGA_DIR_NAME + "all_data_by_genes_transformed.txt";
//        String output = DIR_NAME + "TCGA_BRCA_SampleToMCLModuleToValue_AllSamples_070913.txt";
        String output = DIR_NAME + "TCGA_BRCA_SampleToMCLModuleToValue_AllSamples_HigherFilters_070913.txt";
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(tcgaFileName);
        helper.generateSampleToGeneExpClusters(geneToSampleToValue, 
                                               output, 
                                               modules, 
                                               moduleNames);
    }
    
    @Test
    public void selectMCLModules() throws Exception {
        // For treated samples
        String selectedModules = "6, 3, 4, 8";
        String fileName = DIR_NAME + "CoxPHMCLModules_TreatedSamples_061113.txt";
        String geneListFileName = DIR_NAME + "MCLModulesGenes_TreatedSamples_061113.txt";
        String attributeFileName = DIR_NAME + "MCLModulesAtt_TreatedSamples_061113.txt";
        Map<String, Integer> treatedGeneToModule = selectMCLModules(selectedModules, fileName, geneListFileName, attributeFileName);
        // For all samples
        selectedModules = "33, 2, 46, 37, 34, 15, 28, 31, 39, 45, 42, 29, 5";
        fileName = DIR_NAME + "CoxPHMCLModules_AllSamples_061113.txt";
        geneListFileName = DIR_NAME + "MCLModulesGenes_AllSamples_061113.txt";
        attributeFileName = DIR_NAME + "MCLModulesAtt_AllSamples_061113.txt";
        Map<String, Integer> geneToModule = selectMCLModules(selectedModules, fileName, geneListFileName, attributeFileName);
        Set<String> shared = InteractionUtilities.getShared(treatedGeneToModule.keySet(),
                                                            geneToModule.keySet());
        System.out.println("Genes in treated sample modules: " + treatedGeneToModule.size());
        System.out.println("Genes in all sample modules: " + geneToModule.size());
        System.out.println("Shared genes: " + shared.size());
    }
    
    private Map<String, Integer> selectMCLModules(String moduleText,
                                                  String fileName,
                                                  String geneListFileName,
                                                  String attributeFileName) throws IOException {
        Set<Integer> selectedModules = new HashSet<Integer>();
        if (moduleText != null) {
            String[] modules = moduleText.split(", ");
            for (String module : modules)
                selectedModules.add(new Integer(module));
        }
        Map<String, Integer> geneToModule = new HashMap<String, Integer>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String text = tokens[0].substring("Module".length());
            Integer module = new Integer(text);
            if (selectedModules.size() == 0 || 
                selectedModules.contains(module + 1)) { // Module output in R is 1 based
                String geneText = tokens[tokens.length - 1];
                int index1 = geneText.indexOf("[");
                int index2 = geneText.indexOf("]");
                geneText = geneText.substring(index1 + 1, index2);
                String[] genes = geneText.split(", ");
                for (String gene : genes)
                    geneToModule.put(gene, module);
            }
        }
        fu.close();
        if (geneListFileName != null) {
            fu.setOutput(geneListFileName);
            for (String gene : geneToModule.keySet())
                fu.printLine(gene);
            fu.close();
        }
        if (attributeFileName != null)
            CancerAnalysisUtilitites.outputMapAsCytoscapeAttribute(geneToModule, 
                                                                   attributeFileName,
                                                                   "Module", 
                                                                   Integer.class);
        return geneToModule;
    }
    
    /**
     * Do a MCL based clustering analysis.
     * @throws Exception
     */
    @Test
    public void doMCLClusteringBasedOnPValue() throws Exception {
//        String scoreFileName = DIR_NAME + "CoxPHSurvivalAnalysisBasedOnCNVGenes_061013.txt";
        String scoreFileName = DIR_NAME + "CoxPHSurvivalAnalysisBasedOnCNVGenesLongerDRFS_070313.txt";
        Map<String, Double> geneToPvalue = new HashMap<String, Double>();
        fu.setInput(scoreFileName);
        String line = fu.readLine();
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // 2 for all samples, and 5 for treated samples
            geneToPvalue.put(tokens[0], new Double(tokens[2]));
//            geneToPvalue.put(tokens[0], new Double(tokens[5]));
        }
        fu.close();
        
//        // Add a randomization
//        geneToPvalue = MathUtilities.permutate(geneToPvalue);
        
        String fileName = DIR_NAME + "NTcnMaxVsRefGenesHG19_transformed.txt";
////        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
//        // Use longerFU_DRFS as suggested by John
        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed_longerDRFS.txt";
        // For TCGA output
//        String fileName = TCGA_DIR_NAME + "all_data_by_genes_transformed.txt";
//        String clinFile = DIR_NAME + "TCGA_BRCA.clin.merged.picked.transformed.txt";
        
        TCGAFireHoseOVAnalyzer helper = new TCGAFireHoseOVAnalyzer();
        SurvivalAnalysisResult results = helper.doMCLClusteringBasedOnPValue(geneToPvalue,
                                                                             fileName,
                                                                             clinFile,
                                                                             null,
                                                                             null);
        System.out.println("Error: \n" + results.getError());
        String svResult = CancerAnalysisUtilitites.attachFDRsToSurvivalResults(results.getOutput());
        System.out.println("Result: \n" + svResult);
    }
    
    /**
     * Do a gene-based survival analysis.
     * @throws IOException
     */
    @Test
    public void doGeneBasedSurvivalAnalysis() throws IOException {
        String fileName = DIR_NAME + "NTcnMaxVsRefGenesHG19_transformed.txt";
//        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
        // Use longerFU_DRFS as suggested by John
        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed_longerDRFS.txt";
        
        CancerGeneExpressionCommon loader = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = loader.loadGeneExpAfterAverage(fileName);
        System.out.println("Total genes: " + geneToSampleToValue.size());
        
        // Create a temp score file from sample to gene to value
        String scoreFileName = R3Constants.TEMP_DIR + "treated_score.txt";
        Set<String> treatedPatients = getTreatedPatients();
        loader.generateSampleToGeneToExpValue(geneToSampleToValue, 
                                              treatedPatients,
                                              scoreFileName);
        
//        String scoreFileName = R3Constants.TEMP_DIR + "score.txt";
//        loader.generateSampleToGeneToExpValue(geneToSampleToValue, 
//                                              scoreFileName);
        
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        SurvivalAnalysisResult result = survivalHelper.doSurvivalAnalysis(new File(scoreFileName),
                                                                          new File(clinFile),
                                                                          "coxph",
                                                                          null,
                                                                          "cnv",
                                                                          false);
        System.out.println("Error:\n" + result.getError());
        String svResult = CancerAnalysisUtilitites.attachFDRsToSurvivalResults(result.getOutput());
        System.out.println("Results:\n" + svResult);
    }
    
    /**
     * Check genes contained by gene based CNV file.
     * @throws IOException
     */
    @Test
    public void checkGeneBasedFile() throws IOException {
        String fileName = DIR_NAME + "NTcnMaxVsRefGenesHG19_transformed.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        List<String> geneList = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneList.add(tokens[0]);
        }
        fu.close();
        System.out.println("Total genes in list: " + geneList.size());
        Set<String> geneSet = new HashSet<String>(geneList);
        System.out.println("Total genes in set: " + geneSet.size());
    }
    
    /**
     * Calculate FDRs for a list of p-values
     * @throws IOException
     */
    @Test
    public void calculateFDRs() throws IOException {
        String fileName = DIR_NAME + "pvalues.txt";
        fu.setInput(fileName);
        List<Double> pvalues = new ArrayList<Double>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            Double pvalue = new Double(line);
            pvalues.add(pvalue);
        }
        fu.close();
        Collections.sort(pvalues);
        List<Double> fdrs = MathUtilities.calculateFDRWithBenjaminiHochberg(pvalues);
        System.out.println("pvalue\tFDR");
        for (int i = 0; i < pvalues.size(); i++) {
            Double pvalue = pvalues.get(i);
            Double fdr = fdrs.get(i);
            System.out.println(pvalue + "\t" + fdr);
        }
    }
    
    @Test
    public void discretizeSegmentBasedFile() throws IOException {
        String dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed.txt";
//        String outFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_disc_bug_fix.txt";
        String outFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_0_disc.txt";
        Map<String, Map<String, Double>> sampleToCNVToValue = new CancerGeneExpressionCommon().loadGeneExp(dataFile);
        System.out.println("Total CNV segments: " + sampleToCNVToValue.size());
        Set<String> cnvs = new HashSet<String>();
        for (String sample : sampleToCNVToValue.keySet()) {
            Map<String, Double> cnvToValue = sampleToCNVToValue.get(sample);
            cnvs.addAll(cnvToValue.keySet());
        }
        Map<String, Double> cnvToThreshold = new HashMap<String, Double>();
        for (String cnv : cnvs) {
            List<Double> values = new ArrayList<Double>();
            for (String sample : sampleToCNVToValue.keySet()) {
                Map<String, Double> cnvToValue = sampleToCNVToValue.get(sample);
                values.add(cnvToValue.get(cnv));
            }
            //Double threshold = searchThreshold(values);
            //cnvToThreshold.put(cnv, threshold);
            cnvToThreshold.put(cnv, 0.0d);
        }
        for (String cnv : cnvToThreshold.keySet())
            System.out.println(cnv + "\t" + cnvToThreshold.get(cnv));
        discretize(dataFile, outFile, cnvToThreshold);
    }

    private Double searchThreshold(List<Double> valueList) {
        // Get some summary stats
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (Double v : valueList)
            stat.addValue(v);
        System.out.println("Average: " + stat.getMean());
        System.out.println("Median: " + stat.getPercentile(50.0d));
        Double median = stat.getPercentile(50.0d);
        // Convert value list into set
        Set<Double> valueSet = new HashSet<Double>(valueList);
//        System.out.println("Total different values: " + valueSet.size());
        // Do a sorting
        List<Double> distValueList = new ArrayList<Double>(valueSet);
        Collections.sort(distValueList);
        // Need to find threshold so that the distance between two values closed to this threshold
        // (one up and one down) the the biggest. This is the method used by Alex at CSHL based on
        // email forwarded from Lincoln on June 26, 2013.
        Double maxDiff = Double.MIN_VALUE;
        Double minDistToMedian = Double.MAX_VALUE;
        Double threshold = null;
        for (int i = 0; i < distValueList.size() - 1; i++) {
            Double value1 = distValueList.get(i);
            Double value2 = distValueList.get(i + 1);
            Double diff = value2 - value1;
            Double tmp = (value2 + value1) / 2.0d;
            Double distToMedian = Math.abs(tmp - median);
            // Want to choose a value close to the median
            // so that we can get a even split of samples
            if (diff > maxDiff) {
                threshold = tmp;
                maxDiff = diff;
                minDistToMedian = distToMedian;
            }
            else if (diff == maxDiff && distToMedian < minDistToMedian) { // In case of the same diff, we choose based on distance to median
                threshold = tmp;
                minDistToMedian = distToMedian;
            }
        }
//        System.out.println("Threshold: " + threshold);
        return threshold;
    }
    
    private void discretize(String srcFile,
                            String outFile,
                            Map<String, Double> cnvToTreshold) throws IOException {
        fu.setInput(srcFile);
        fu.setOutput(outFile);
        String line = fu.readLine();
        String[] cnvs = line.split("\t");
        fu.printLine(line); // Print out the header
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            builder.append(tokens[0]);
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("NA"))
                    builder.append("\t").append(tokens[i]);
                else {
                    Double value = new Double(tokens[i]);
                    Double threshold = cnvToTreshold.get(cnvs[i]);
                    if (value <= threshold)
                        builder.append("\t0");
                    else
                        builder.append("\t1");
                }
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * Do survival analysis based on CNV segments.
     * @throws Exception
     */
    @Test
    public void doSegmentBasedSurvivalAnalysis() throws Exception {
        String dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed.txt";
        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed_longerDRFS.txt";
//        String dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_disc.txt";
//        String dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_disc_bug_fix.txt";
//        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
        
        SurvivalAnalysisHelper helper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        SurvivalAnalysisResult result = helper.doSurvivalAnalysis(new File(dataFile), 
                                                                  new File(clinFile), 
                                                                  "coxph", 
                                                                  null, 
                                                                  "BreastAnt", 
                                                                  false);
        System.out.println("CoxPH survival analysis using all samples:");
        System.out.println("Error: \n" + result.getError());
        System.out.println("Result:\n" + result.getOutput());
        
        System.out.println("\nCoxPH survival analysis using treated samples:");
        dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_treated.txt";
        result = helper.doSurvivalAnalysis(new File(dataFile), 
                                           new File(clinFile), 
                                           "coxph", 
                                           null, 
                                           "BreastAnt", 
                                           false);
        System.out.println("Error: \n" + result.getError());
        System.out.println("Result:\n" + result.getOutput());
        
        // Check some significant fragments
        // However, the following code cannot work since kaplan-meier
        // is used to pick up the mean as the pivotal value to split
        // samples. This cannot work for CNV data set.
        // Move the analysis into R directly.
//        result = helper.doSurvivalAnalysis(new File(dataFile), 
//                                           new File(clinFile), 
//                                           "kaplan-meier", 
//                                           "143", 
//                                           "A194", 
//                                           false);
//        System.out.println("Error: \n" + result.getError());
//        System.out.println("Result:\n" + result.getOutput());
    }
    
    /**
     * Process NT20k file so that it can be loaded into R for hierarchical cluster analysis.
     * The file is to be used by HeatMapDrawer.R.
     * @throws Exception
     */
    @Test
    public void processNT20kArchNormFile() throws Exception {
        String inFile = DIR_NAME + "NT/NT20kArchNormIncidenceQCd.txt";
        String outFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed.txt";
        fu.setInput(inFile);
        fu.setOutput(outFile);
        String line = fu.readLine();
        fu.printLine("Sample\t" + line);
        // Remove three lines
        for (int i = 0; i < 4; i++)
            fu.readLine();
        while ((line = fu.readLine()) != null) {
            fu.printLine(line);
        }
        fu.close();
        // Output treated samples only
        Set<String> treatedPatients = getTreatedPatients();
        String outFile1 = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed_treated.txt";
        fu.setOutput(outFile1);
        fu.setInput(outFile);
        line = fu.readLine();
        fu.printLine(line);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (treatedPatients.contains(tokens[0]))
                fu.printLine(line);
        }
        fu.close();
    }

    private Set<String> getTreatedPatients() throws IOException {
        String line;
        // Created a file contains anthracycline treated patients only
        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
        Set<String> treatedPatients = new HashSet<String>();
        fu.setInput(clinFile);
        line = fu.readLine();
        String[] tokens = line.split("\t");
        // Get the index wanted
        int wantedIndex = 0;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("TREAT")) {
                wantedIndex = i;
                break;
            }
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String sample = tokens[0];
            String treat = tokens[wantedIndex];
            if (treat.equals("1"))
                treatedPatients.add(sample);
        }
        fu.close();
        return treatedPatients;
    }
    
    @Test
    public void checkSamples() throws Exception {
        String clinFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
        Set<String> allClinSamples = getSamples(clinFile);
        System.out.println("Total samples in clin: " + allClinSamples.size());
        Set<String> treatedSamples = getTreatedPatients();
        System.out.println("Treated samples: " + treatedSamples.size());
        
        // Samples in data files
        String dataFile = DIR_NAME + "NT20kArchNormIncidenceQCd_transformed.txt";
        Set<String> segmentSamples = getSamples(dataFile);
        System.out.println("Total samples in segments: " + segmentSamples.size());
        segmentSamples.retainAll(allClinSamples);
        System.out.println("Having clin: " + segmentSamples.size());
        segmentSamples.retainAll(treatedSamples);
        System.out.println("Treated: " + segmentSamples.size());
    }
    
    private Set<String> getSamples(String file) throws IOException {
        Set<String> samples = new HashSet<String>();
        fu.setInput(file);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            samples.add(tokens[0]);
        }
        fu.close();
        return samples;
    }
    
    @Test
    public void checkNT20kArchNormFiles() throws Exception {
        String fileName = DIR_NAME + "NT/NT20kArchNormIncidenceQCd.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        List<String> segments = new ArrayList<String>();
        for (String token : tokens)
            segments.add(token);
        System.out.println("Total segments: " + tokens.length);
        fu.close();
        
        fileName = DIR_NAME + "NT20kArchNormIncidence.txt";
        fu.setInput(fileName);
        line = fu.readLine();
        tokens = line.split("\t");
        List<String> segments1 = new ArrayList<String>();
        for (String token : tokens)
            segments1.add(token);
        System.out.println("Total segments: " + segments1.size());
        segments1.removeAll(segments);
        System.out.println("Extra segments: " + segments1.size());
        for (String seg : segments1)
            System.out.println(seg);
        fu.close();
    }
    
    /**
     * Make a little change to the clinical file so that NT ids are the same in other CNV files.
     * @throws IOException
     */
    @Test
    public void processClinFile() throws IOException {
        String inFile = DIR_NAME + "McNEAT_Clin_Info.txt";
        String outFile = DIR_NAME + "McNEAT_Clin_Info_transformed.txt";
        fu.setInput(inFile);
        fu.setOutput(outFile);
        String line = fu.readLine();
        int index = line.indexOf("\t");
        fu.printLine(line.substring(index + 1));
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[1];
            while (sample.length() < 3)
                sample = "0" + sample;
            tokens[1] = "NT" + sample;
            // But we don't want the first column
            List<String> list = new ArrayList<String>();
            for (int i = 1; i < tokens.length; i++)
                list.add(tokens[i]);
            fu.printLine(StringUtils.join(list.iterator(), "\t"));
        }
        fu.close();
    }
}
