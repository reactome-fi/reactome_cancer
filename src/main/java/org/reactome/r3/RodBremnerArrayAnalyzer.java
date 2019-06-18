/*
 * Created on Mar 10, 2009
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.StatisticalSummary;
import org.apache.commons.math.stat.inference.TestUtils;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;

/**
 * This class is used to do array related analysis for Rob Bremner.
 * @author wgm
 *
 */
public class RodBremnerArrayAnalyzer {
    private final String DATA_SET = "datasets/Rod_IFNG/arrays/";
    private final String SOURCE_FILE = DATA_SET + "GEX Data/TXT tables from BeadStudio/8grps-24smpls_norm_bkgd/Sample_Gene_Norm_Cubic_Profille.txt";
    private final String RESULT_DIR = DATA_SET + "IFNGArrayResults/";
    private FileUtility fu = new FileUtility();
    
    public RodBremnerArrayAnalyzer() {
    }
    
    @Test
    public void checkDiffGenesValidationForMohamed() throws Exception {
        String dirName = DATA_SET + "mohamed/";
        String fileName = dirName + "RT_qPCRValidationsPValue.txt";
        Map<String, Double> geneToPvalue = new HashMap<String, Double>();
        int index = 3;
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            geneToPvalue.put(tokens[4], new Double(tokens[index]));
        }
        fu.close();
        List<String> validatedGenes = new ArrayList<String>();
        List<String> invalidatedGenes = new ArrayList<String>();
        for (String gene : geneToPvalue.keySet()) {
            Double pvalue = geneToPvalue.get(gene);
            if (pvalue.equals(Double.NaN))
                invalidatedGenes.add(gene);
            else if (pvalue <= 0.05)
                validatedGenes.add(gene);
            else
                invalidatedGenes.add(gene);
        }
        fileName = dirName + "BRG1_IFN_vs_SW13_up_genes_BH_0619_cut_02.txt";
        List<String> diffUpGenes = loadGenes(fileName);
        // Check how many genes in the diffUpGenes
        System.out.println("Total validated genes: " + validatedGenes.size());
        List<String> copy = new ArrayList<String>(validatedGenes);
        validatedGenes.retainAll(diffUpGenes);
        System.out.println("validated genes in diff up list: " + validatedGenes.size());
        // Want to find what genes are not in our diff list
        copy.removeAll(validatedGenes);
        System.out.println("validated genes not in the diff list: " + copy.size() + " " + copy);
        System.out.println("Total invalidated genes: " + invalidatedGenes.size());
        invalidatedGenes.retainAll(diffUpGenes);
        System.out.println("Invalidated genes in diff up list: " + invalidatedGenes.size());
    }
    
    @Test
    public void checkRTqPCRValidationResults() throws Exception {
        String fileName = DATA_SET + "mohamed/RT_qPCRValidations.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        line = fu.readLine();
        // For t-value calculate
        ValidationResult r1 = new ValidationResult();
        ValidationResult r2 = new ValidationResult();
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            // As a reference
            r1.mean = new Double(tokens[1]);
            r1.sd = new Double(tokens[6]);
            builder.append(gene);
            for (int i = 2; i < 5; i++) {
                r2.mean = new Double(tokens[i]);
                r2.sd = new Double(tokens[i + 5]);
                double pvalue = TestUtils.tTest(r1, r2) / 2.0;
                if (r2.mean < r1.mean)
                    pvalue = 1.0 - pvalue;
                builder.append("\t").append(pvalue);
            }
            System.out.println(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    @Test
    public void checkSharedGenes() throws Exception {
        String dirName = DATA_SET + "mohamed/";
        String fileName = dirName + "BRG1_vs_SW13_up_genes_BH_0619_cut_03.txt";
        List<String> brg1UpGenes = loadGenes(fileName);
        fileName = dirName + "BRG1_vs_SW13_down_genes_BH_0619_cut_03.txt";
        List<String> brg1DownGenes = loadGenes(fileName);
        fileName = dirName + "BRG1_IFN_vs_SW13_up_genes_BH_0619_cut_03.txt";
        List<String> brg1IfnUpGenes = loadGenes(fileName);
        fileName = dirName + "BRG1_IFN_vs_SW13_down_genes_BH_0619_cut_03.txt";
        List<String> brg1IfnDownGenes = loadGenes(fileName);
        fileName = dirName + "IFN_vs_SW13_up_genes_BH_0619_cut_03.txt";
        List<String> ifnUpGenes = loadGenes(fileName);
        fileName = dirName + "IFN_vs_SW13_down_genes_BH_0619_cut_03.txt";
        List<String> ifnDownGenes = loadGenes(fileName);
        // For up genes
        System.out.println("Up genes: ");
        checkSharedGenes(brg1UpGenes, brg1IfnUpGenes, ifnUpGenes);
        System.out.println("\nDown genes:");
        checkSharedGenes(brg1DownGenes, brg1IfnDownGenes, ifnDownGenes);
    }
    
    @Test
    public void countShareInPairwiseGenesForDax() throws IOException {
        String input = DATA_SET + "dax/CountsForCategoriesOfGenes_Cutoff_Dax_cut_05.txt";
        String output = DATA_SET + "dax/CountsForCategoriesOfGenes_Cutoff_Dax_cut_05_withcounts.txt";
        fu.setInput(input);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(output);
        String line = fu.readLine();
        // Title line
        StringBuilder builder = new StringBuilder();
        builder.append(line);
        builder.append("\tBRG1hi\tBRG1lo\tSTATlo\tBRG1hi+BRG1lo\tBRG1hi+BRG1lo+STAT10");
        outFu.printLine(builder.toString());
        int brgHi = 0;
        int brgLo = 0;
        int statLo = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // BRGhi
            brgHi = countShares(tokens, 1, 5);
            brgLo = countShares(tokens, 5, 13);
            statLo = countShares(tokens, 13, 17);
            builder.setLength(0);
            builder.append(line).append("\t");
            outputShares(builder, brgHi);
            builder.append("\t");
            outputShares(builder, brgLo);
            builder.append("\t");
            outputShares(builder, statLo);
            builder.append("\t");
            // Check combination: need to flip the sign of brgLow since we want to see the sharing
            // between up genes in brgHi and down genes in brgLo
            int shared = countShares(brgHi, -brgLo);
            outputShares(builder, shared);
            builder.append("\t");
            shared = countShares(brgHi, -brgLo, -statLo);
            outputShares(builder, shared);
            outFu.printLine(builder.toString());
        }
        outFu.close();
        fu.close();
    }
    
    private int countShares(String[] tokens, 
                            int startIndex, 
                            int endIndex) {
        Set<String> set = new HashSet<String>();
        int count = 0;
        for (int i = startIndex; i < endIndex; i++) {
            set.add(tokens[i]);
            count += Integer.parseInt(tokens[i]);
        }
        if (set.contains("1") && set.contains("-1"))
            return Integer.MIN_VALUE; // Used as a mark
        return count;
    }
    
    private int countShares(int ... counts) {
        // Make sure only one case is possible
        Set<String> cases = new HashSet<String>();
        int total = 0;
        for (int count : counts) {
            if (count == Integer.MIN_VALUE)
                return Integer.MIN_VALUE;
            if (count > 0)
                cases.add("+");
            else if (count < 0)
                cases.add("-");
            total += count;
        }
        if (cases.size() == 2)
            return Integer.MIN_VALUE;
        return total;
    }
    
    private void outputShares(StringBuilder builder,
                              int counter) {
        if (counter == Integer.MIN_VALUE)
            builder.append("NA");
        else
            builder.append(counter);
    }
                              
    
    /**
     * Check overlapping for up/down genes based on pair-wise analysis.
     * @throws Exception
     */
    @Test
    public void checkShareInPairwiseGenesForDax() throws Exception {
        String[] fileNames = new String[] {
                "2FBRG1_vs_2FD3_down_genes_RAWP_0210_cut_05.txt",          "2FD3_siBRG_vs_2FGFP_down_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FD3_siCONT_down_genes_RAWP_0210_cut_05.txt",       "2FD3_siBRG_vs_2FGFP_siCONT_down_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FD3_siCONT_up_genes_RAWP_0210_cut_05.txt",     "2FD3_siBRG_vs_2FGFP_siCONT_up_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FD3_up_genes_RAWP_0210_cut_05.txt",            "2FD3_siBRG_vs_2FGFP_up_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FGFP_down_genes_RAWP_0210_cut_05.txt",         "U3A_vs_2FD3_down_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FGFP_siCONT_down_genes_RAWP_0210_cut_05.txt",      "U3A_vs_2FD3_siCONT_down_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FGFP_siCONT_up_genes_RAWP_0210_cut_05.txt",        "U3A_vs_2FD3_siCONT_up_genes_RAWP_0210_cut_05.txt",
                "2FBRG1_vs_2FGFP_up_genes_RAWP_0210_cut_05.txt",           "U3A_vs_2FD3_up_genes_RAWP_0210_cut_05.txt",
                "2FD3_siBRG_vs_2FD3_down_genes_RAWP_0210_cut_05.txt",      "U3A_vs_2FGFP_down_genes_RAWP_0210_cut_05.txt",
                "2FD3_siBRG_vs_2FD3_siCONT_down_genes_RAWP_0210_cut_05.txt",   "U3A_vs_2FGFP_siCONT_down_genes_RAWP_0210_cut_05.txt",
                "2FD3_siBRG_vs_2FD3_siCONT_up_genes_RAWP_0210_cut_05.txt",     "U3A_vs_2FGFP_siCONT_up_genes_RAWP_0210_cut_05.txt",
                "2FD3_siBRG_vs_2FD3_up_genes_RAWP_0210_cut_05.txt",        "U3A_vs_2FGFP_up_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FGFP_siCONT_up_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FGFP_siCONT_down_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FGFP_up_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FGFP_down_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FD3_siCONT_up_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FD3_siCONT_down_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FD3_up_genes_RAWP_0210_cut_05.txt",
                "2FGFP_siBRG_vs_2FD3_down_genes_RAWP_0210_cut_05.txt",
        };
        String dirName = DATA_SET + "dax/";
        // Get all genes used
        Set<String> allGenes = fu.loadInteractions(dirName + "Gene_List_021010.txt");
        List<String> allGeneList = new ArrayList<String>(allGenes);
        Collections.sort(allGeneList);
        Map<String, Set<String>> fileToGenes = new HashMap<String, Set<String>>();
        for (String file : fileNames) {
            Set<String> genes = fu.loadInteractions(dirName + file);
            fileToGenes.put(file, genes);
        }
        System.out.println("Total files: " + fileToGenes.size());
//        // Check BRG1 up genes
//        for (String file1 : fileToGenes.keySet()) {
//            if (!file1.startsWith("2FBRG1"))
//                continue;
//            Set<String> genes1 = fileToGenes.get(file1);
//            for (String file2 : fileToGenes.keySet()) {
//                if (!file2.startsWith("2FBRG1"))
//                    continue;
//                if (file1.equals(file2))
//                    continue;
//                // Want to compare up vs up, down vs down
//                if (file1.contains("up_genes") && file2.contains("down_genes"))
//                    continue;
//                if (file1.contains("down_genes") && file2.contains("up_genes"))
//                    continue;
//                Set<String> genes2 = fileToGenes.get(file2);
//                Set<String> shared = InteractionUtilities.getShared(genes1, genes2);
//                Double pvalue = MathUtilities.calculateHypergeometricPValue(7811, // Total valid genes 
//                                                                            genes1.size(), 
//                                                                            genes2.size(), 
//                                                                            shared.size());
//                System.out.println(file1 + "\t" + file2 + "\t" + pvalue);
//            }
//        }
        StringBuilder builder = new StringBuilder();
        String fileName = dirName + "CategoriesOfGenes_Cutoff_Dax_05.txt";
        String inducedFileName = dirName + "CategoriesOfGenes_Cutoff_Dax_05_induced.txt";
        String repressedFileName = dirName + "CategoriesOfGenes_Cutoff_Dax_05_repressed.txt";
        fu.setOutput(fileName);
        builder.append("GeneName\t2fBRG v 2fD3\t2fBRG v siCONT_2fD3\t2fBRG v 2fGFP\t2fBRG v siCONT_2fGFP\t" + 
                       "siBRG_2fD3 v 2fD3\tsiBRG_2fD3 v siCONT_2fD3\tsiBRG_2fD3 v 2fGFP\tsiBRG_2fD3 v siCONT_2fGFP\t" +
                       "siBRG_2fGFP v 2fD3\tsiBRG_2fGFP v siCONT_2fD3\tsiBRG_2fGFP v 2fGFP\tsiBRG_2fGFP v siCONT_2fGFP\t" + 
                       "U3A v 2fD3\tU3A v siCONT_2fD3\tU3A v 2fGFP\tU3A v siCONT_2fGFP\t" + 
                       "Up 4/4 BRG1hi\tDown 8/8 BRG1lo\tUp 4/4 BRG1hi & Down 8/8 BRG1lo\t" +
                       "Up 3/4 BRG1hi\tDown 7/8 BRG1lo\tUp 3/4 BRG1hi & Down 7/8 BRG1lo\t" + 
                       "Down in 4/4 STAT1lo\tDown in 3/4 STAT1lo\t" + 
                       "Up 4/4 BRG1hi, down 4/4 STAT1lo\tDown 8/8 BRG1lo, down 4/4 STAT1lo\tUp 4/4 BRG1hi & Down 8/8 BRG1lo, & down 4/4 STAT1lo\t" +
                       "Up 3/4 BRG1hi, down 3/4 STAT1lo\tDown 7/8 BRG1lo, down 3/4 STAT1lo\tUp 3/4 BRG1hi & Down 7/8 BRG1lo, & down 3/4 STAT1lo\t" +
                       "Down 4/4 BRG1 hi\tUp 8/8 BRG1lo\tDown 4/4 BRG1hi & Up 8/8 BRG1lo\t" + 
                       "Down 3/4 BRG1 hi\tUp 7/8 BRG1lo\tDown 3/4 BRG1hi & Up 7/8 BRG1lo");
        fu.printLine(builder.toString());
        // Used to count the types
        List<String> list = new ArrayList<String>();
        String[] samplePairs = new String[] {
                // Check BRG1 and controls (up genes)
                "2FBRG1_vs_2FD3",
                "2FBRG1_vs_2FD3_siCONT",
                "2FBRG1_vs_2FGFP",
                "2FBRG1_vs_2FGFP_siCONT",
                // Check siBRG and controls (down genes)
                "2FD3_siBRG_vs_2FD3",
                "2FD3_siBRG_vs_2FD3_siCONT",
                "2FD3_siBRG_vs_2FGFP",
                "2FD3_siBRG_vs_2FGFP_siCONT",
                "2FGFP_siBRG_vs_2FD3",
                "2FGFP_siBRG_vs_2FD3_siCONT",
                "2FGFP_siBRG_vs_2FGFP",
                "2FGFP_siBRG_vs_2FGFP_siCONT",
                // U3A and controls (down genes)
                "U3A_vs_2FD3",
                "U3A_vs_2FD3_siCONT",
                "U3A_vs_2FGFP",
                "U3A_vs_2FGFP_siCONT"
        };
        for (String gene : allGeneList) {
            gene = gene.trim();
            builder.setLength(0);
            builder.append(gene);
            list.clear();
            for (String samplePair : samplePairs) {
                getSharingResult(gene, fileToGenes, samplePair, list);
            }
            for (String share : list)
                builder.append("\t").append(share);
            // Check sharing for stringent BRG1-induced genes
            int brgUpInHi = countSharingResult(list, 0, 3, "+1");
            if (brgUpInHi == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            int brgDownInLo = countSharingResult(list, 4, 11, "-1");
            if (brgDownInLo == 8)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgUpInHi == 4 && brgDownInLo == 8)
                builder.append("\t1");
            else
                builder.append("\t0");
            // Less stringent
            if (brgUpInHi >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgDownInLo >= 7)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgUpInHi >= 3 && brgDownInLo >= 7)
                builder.append("\t1");
            else
                builder.append("\t0");
            int downStat1Lo = countSharingResult(list, 12, 15, "-1");
            if (downStat1Lo == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (downStat1Lo >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            // Stringent BRG1 and STAT1
            if (brgUpInHi == 4 && downStat1Lo == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgDownInLo == 8 && downStat1Lo == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgUpInHi == 4 && brgDownInLo == 8 && downStat1Lo == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            // Less stringent BRG1 and STAT1
            if (brgUpInHi >= 3 && downStat1Lo >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgDownInLo >= 7 && downStat1Lo >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (brgUpInHi >= 3 && brgDownInLo >= 7 && downStat1Lo >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            // Stringent BRG1 repressed genes
            int downBRGInHi = countSharingResult(list, 0, 3, "-");
            int upBRGInLo = countSharingResult(list, 4, 11, "+");
            if (downBRGInHi == 4)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (upBRGInLo == 8)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (downBRGInHi == 4 && upBRGInLo == 8)
                builder.append("\t1");
            else
                builder.append("\t0");
            // Less stringent BRG1 repressed genes
            if (downBRGInHi >= 3)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (upBRGInLo >= 7)
                builder.append("\t1");
            else
                builder.append("\t0");
            if (downBRGInHi >= 3 && upBRGInLo >= 7)
                builder.append("\t1");
            else
                builder.append("\t0");
            fu.printLine(builder.toString());
        }
        fu.close();
        // Want to create a pair-wise overlapping analsis to check the number and significants of shared genes
        // BRG1 induced genes
        checkSharedGenes(allGeneList, 
                         fileToGenes, 
                         inducedFileName,
                         samplePairs,
                         true);
        // BRG1 repressed genes
        checkSharedGenes(allGeneList, 
                         fileToGenes, 
                         repressedFileName,
                         samplePairs,
                         false);
    }

    private void checkSharedGenes(List<String> allGeneList,
                                  Map<String, Set<String>> fileToGenes,
                                  String inducedFileName,
                                  String[] samplePairs,
                                  boolean isForInduced) throws IOException, MathException {
        fu.setOutput(inducedFileName);
        StringBuilder builder = new StringBuilder();
        for (String samplePair : samplePairs) {
            builder.append("\t").append(samplePair);
        }
        fu.printLine(builder.toString());
        boolean isUp;
        for (int i = 0; i < samplePairs.length; i++) {
            builder.setLength(0);
            String sample1 = samplePairs[i];
            builder.append(sample1);
            if (isForInduced) {
                isUp = (i < 4) ? true : false;
            }
            else {
                isUp = (i < 4) ? false : true;
            }
            Set<String> genes1 = getDiffGeneSets(sample1,
                                                 isUp,
                                                 fileToGenes);
            for (int j = 0; j < samplePairs.length; j++) {
                if (i == j) {
                    builder.append("\t");
                    continue;
                }
                if (isForInduced) {
                    isUp = (j < 4) ? true : false;
                }
                else {
                    isUp = (j < 4) ? false : true;
                }
                String sample2 = samplePairs[j];
                Set<String> genes2 = getDiffGeneSets(sample2,
                                                     isUp,
                                                     fileToGenes);
                Set<String> shared = InteractionUtilities.getShared(genes1, genes2);
                double pvalue = MathUtilities.calculateHypergeometricPValue(allGeneList.size(),
                                                                            genes1.size(),
                                                                            genes2.size(),
                                                                            shared.size());
                builder.append("\t").append(shared.size()).append(" (").append(String.format("%2.2e", pvalue)).append(")");
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    private int countSharingResult(List<String> list,
                                   int startIndex,
                                   int endIndex,
                                   String targetSymbol) {
        int count = 0;
        for (int i = startIndex; i <= endIndex; i++) {
            String tmp = list.get(i);
            if (tmp.equals(targetSymbol))
                count++;
        }
        return count;
    }
    
    private void getSharingResult(String gene,
                                  Map<String, Set<String>> fileToGenes,
                                  String samplePair,
                                  List<String> list) {
        Set<String> downGenes = getDiffGeneSets(samplePair, 
                                                false,
                                                fileToGenes);
        Set<String> upGenes = getDiffGeneSets(samplePair,
                                              true,
                                              fileToGenes);
        if (downGenes.contains(gene)) 
            list.add("-1");
        else if (upGenes.contains(gene)) 
            list.add("+1");
        else
           list.add("0");
    }
                                      
    
    private Set<String> getDiffGeneSets(String samplePair,
                                        boolean isForUp,
                                        Map<String, Set<String>> fileToGenes) {
        String key = samplePair + "_" + (isForUp ? "up_" : "down_") + "genes";
        for (String file : fileToGenes.keySet()) {
            if (file.startsWith(key))
                return fileToGenes.get(file);
        }
        System.err.println(samplePair + " is not correct!");
        return null;
    }
    
    @Test
    public void checkBRG1GenesInTwoDatasets() throws Exception {
        int totalGenes = 14238; // Use the samllest number
        String daxDir = DATA_SET + "dax/";
        String fileName = daxDir + "BRG_vs_SiCONT_GFP_up_genes_BH.txt";
        String list1Name = "Dax BRG1";
        List<String> list1 = loadGenes(fileName);
        String mohamedDir = DATA_SET + "mohamed/";
        fileName = mohamedDir + "BRG1_vs_SW13_up_genes_BH.txt";
        String list2Name = "Mohamed BRG1";
        List<String> list2 = loadGenes(fileName);
        System.out.println("\nOverlapping for BRG1 up between Dax and Mohamed:");
        checkSharedGenes(list1, 
                         list1Name,
                         list2,
                         list2Name,
                         totalGenes);
        fileName = daxDir + "BRG_vs_SiCONT_GFP_down_genes_BH.txt";
        list1Name = "Dax BRG1";
        list1 = loadGenes(fileName);
        fileName = mohamedDir + "BRG1_vs_SW13_down_genes_BH.txt";
        list2Name = "Mohamed BRG1";
        list2 = loadGenes(fileName);
        System.out.println("\nOverlapping for BRG1 down between Dax and Mohamed:");
        checkSharedGenes(list1, 
                         list1Name,
                         list2,
                         list2Name,
                         totalGenes);
    }
    
    @Test
    public void checkSharedGenesInDaxDataset() throws Exception {
        int totalGenes = 16971;
        String dirName = DATA_SET + "dax/";
//        String fileName = dirName + "BRG1_up_genes.txt";
//        List<String> brg1UpGenes = loadGenes(fileName);
//        fileName = DATA_SET + "mohamed/BRG1_up_genes.txt";
//        List<String> mohamedBrg1UpGenes = loadGenes(fileName);
//        List<String> shared = getShared(brg1UpGenes, mohamedBrg1UpGenes);
//        double pvalue = calculatePValueForSharing(totalGenes, 
//                                                  brg1UpGenes,
//                                                  mohamedBrg1UpGenes, 
//                                                  shared);
//        System.out.println("BRG1 up genes in Dax: " + brg1UpGenes.size());
//        System.out.println("BRG1 up genes in Mohamed: " + mohamedBrg1UpGenes.size());
//        System.out.println("Shared BRG1 up genes between Dax and Mohamed: " + shared.size() + 
//                           " (" + pvalue + ")");
        // Compare BRG1, stat1 and siBRG
        String fileName = dirName + "BRG_vs_SiCONT_GFP_up_genes_BH_0619_cut_03.txt";
        String list1Name = "BRG1";
        List<String> list1 = loadGenes(fileName);
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_down_genes_BH_0619_cut_03.txt";
        String list2Name = "siBRG1";
        List<String> list2 = loadGenes(fileName);
        fileName = dirName + "stat1_vs_SiCONT_2FD3_down_genes_BH_0619_cut_03.txt";
        String list3Name = "stat1";
        List<String> list3 = loadGenes(fileName);
        System.out.println("\nOverlapping for BRG1 up and stat1 down genes:");
        checkSharedGenes(list1, 
                         list1Name,
                         list3, 
                         list3Name,
                         totalGenes);
        System.out.println("\nOverlapping for siBRG1 down and stat1 down genes:");
        checkSharedGenes(list2, 
                         list1Name,
                         list3, 
                         list3Name,
                         totalGenes);
        // All three sharing
        List<String> all = getShared(list1, list2, list3);
        System.out.println("\nShared in BRG1 up, stat1 down and siBRG1 down: " + all.size());
        if (all.size() < 50)
            System.out.println(getStringForCollection(all));
        fileName = dirName + "BRG_vs_SiCONT_GFP_down_genes_BH_0619_cut_03.txt";
        list1 = loadGenes(fileName);
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_up_genes_BH_0619_cut_03.txt";
        list2 = loadGenes(fileName);
        fileName = dirName + "stat1_vs_SiCONT_2FD3_up_genes_BH_0619_cut_03.txt";
        list3 = loadGenes(fileName);
        System.out.println("\nOverlapping for BRG1 down and stat1 up genes:");
        checkSharedGenes(list1, 
                         list1Name,
                         list3, 
                         list3Name,
                         totalGenes);
        System.out.println("\nOverlapping for siBRG1 up and stat1 up genes:");
        checkSharedGenes(list2, 
                         list1Name,
                         list3, 
                         list3Name,
                         totalGenes);
        // All three sharing
        all = getShared(list1, list2, list3);
        System.out.println("\nShared in BRG1 down, stat1 up and siBRG1 up: " + all.size());
        if (all.size() < 50)
            System.out.println(getStringForCollection(all));
        // Compare siBRG in two conditions
        // down
        fileName = dirName + "siBRG_vs_SiCONT_2FD3_down_genes_BH_0619_cut_03.txt";
        list1 = loadGenes(fileName);
        list1Name = "siBRG";
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_down_genes_BH_0619_cut_03.txt";
        list2 = loadGenes(fileName);
        list2Name = "GFP-siBRG";
        System.out.println("\nOverlapping for siBRG1 down genes:");
        checkSharedGenes(list1, 
                         list1Name, 
                         list2, 
                         list2Name, 
                         totalGenes);
        // Up
        fileName = dirName + "siBRG_vs_SiCONT_2FD3_up_genes_BH_0619_cut_03.txt";
        list1 = loadGenes(fileName);
        list1Name = "siBRG";
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_up_genes_BH_0619_cut_03.txt";
        list2 = loadGenes(fileName);
        list2Name = "GFP-siBRG";
        System.out.println("\nOverlapping for siBRG1 up genes:");
        checkSharedGenes(list1,
                         list1Name, 
                         list2, 
                         list2Name,
                         totalGenes);
        // Compare siBRG vs BRG
        // BRG1 up and siBRG down 
        fileName = dirName + "BRG_vs_SiCONT_GFP_up_genes_BH_0619_cut_03.txt";
        list1 = loadGenes(fileName);
        list1Name = "BRG1";
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_down_genes_BH_0619_cut_03.txt";
        list2 = loadGenes(fileName);
        list2Name = "siBRG";
        System.out.println("\nOveraling for BRG1 up and siBRG down genes:");
        checkSharedGenes(list1,
                         list1Name, 
                         list2, 
                         list2Name,
                         totalGenes);
        // BRG1 down and siBRG up
        fileName = dirName + "BRG_vs_SiCONT_GFP_down_genes_BH_0619_cut_03.txt";
        list1 = loadGenes(fileName);
        list1Name = "BRG1";
        fileName = dirName + "GFP_siBRG_vs_GFP_SiCONT_GFP_up_genes_BH_0619_cut_03.txt";
        list2 = loadGenes(fileName);
        list2Name = "siBRG";
        System.out.println("\nOveraling for BRG1 down and siBRG up genes:");
        checkSharedGenes(list1,
                         list1Name, 
                         list2, 
                         list2Name,
                         totalGenes);
    }
    
    private void checkSharedGenes(List<String> list1,
                                  String list1Name,
                                  List<String> list2,
                                  String list2Name,
                                  int totalGene) throws MathException {
        List<String> shared = getShared(list1, list2);
        double pvalue = calculatePValueForSharing(totalGene, 
                                                  list1, 
                                                  list2, 
                                                  shared);
        System.out.println(list1Name + ": " + list1.size());
        System.out.println(list2Name + ": " + list2.size());
        System.out.println("shared: " + shared.size() + " (pvalue: " + pvalue + ")");
        if (shared.size() < 50)
            System.out.println(getStringForCollection(shared));
    }

    private void checkSharedGenes(List<String> brg1UpGenes,
                                  List<String> brg1IfnUpGenes,
                                  List<String> ifnUpGenes) throws MathException {
        int totalGenes = 14238;
        //int totalGenes = 7812;
        System.out.println("BRG1 genes: " + brg1UpGenes.size());
        System.out.println("IFN1 genes: " + ifnUpGenes.size());
        System.out.println("BRG1/IFN genes: " + brg1IfnUpGenes.size());       
        // Shared between BRG1 and IFN1
        List<String> shared = getShared(brg1UpGenes, ifnUpGenes);
        // Calculate p-value
        double pvalue = calculatePValueForSharing(totalGenes, 
                                                  brg1UpGenes,
                                                  ifnUpGenes, 
                                                  shared);
        System.out.println("Shared between BRG1 and IFN1: " + shared.size() + " (pvalue: " + pvalue + ")");
        if (shared.size() < 50)
            System.out.println(getStringForCollection(shared));
        shared = getShared(brg1UpGenes, brg1IfnUpGenes);
        pvalue = calculatePValueForSharing(totalGenes,
                                           brg1UpGenes, 
                                           brg1IfnUpGenes,
                                           shared); 
        if (pvalue < 0.0)
            pvalue = 0.0;
        System.out.println("Shared between BRG1 and BRG1/IFN: " + shared.size() + " (pvalue: " + pvalue + ")");
        if (shared.size() < 50)
            System.out.println(getStringForCollection(shared));
        shared = getShared(ifnUpGenes, brg1IfnUpGenes);
        pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes,
                                                             ifnUpGenes.size(),
                                                             brg1IfnUpGenes.size(), 
                                                             shared.size());
        if (pvalue < 0.0)
            pvalue = 0.0;
        System.out.println("Shared between IFN and BRG1/IFN: " + shared.size() + " (pvalue: " + pvalue + ")");
        if (shared.size() < 50)
            System.out.println(getStringForCollection(shared));
        shared = getShared(shared, brg1UpGenes);
        System.out.println("Shared among BRG1, IFN, and BRG1/IFN: " + shared.size());
        if (shared.size() < 50)
            System.out.println(getStringForCollection(shared));
    }

    private double calculatePValueForSharing(int totalGenes,
                                             List<String> brg1UpGenes,
                                             List<String> ifnUpGenes, 
                                             List<String> shared) throws MathException {
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalGenes, 
                                                                    brg1UpGenes.size(),
                                                                    ifnUpGenes.size(), 
                                                                    shared.size());
        if (pvalue < 0.0)
            pvalue = 0.0;
        return pvalue;
    }
    
    private String getStringForCollection(Collection<String> c) {
        String rtn = c.toString();
        return rtn.substring(1, rtn.length() - 1);
    }
    
    private List<String> getShared(List<String>... lists) {
        List<String> copy = new ArrayList<String>(lists[0]);
        for (int i = 1; i < lists.length; i++) 
            copy.retainAll(lists[i]);
        return copy;
    }
    
    private List<String> loadGenes(String fileName) throws IOException {
        List<String> list = new ArrayList<String>();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null)
            list.add(line);
        fu.close();
        return list;
    }
    
    /**
     * This method is used to load a listed genes generated from BioConductor.
     * @return
     * @throws IOException
     */
    public List<String> loadDiffGenesFromFile() throws IOException {
        String fileName = RESULT_DIR + "DiffGenes.txt";
        fu.setInput(fileName);
        String line = null;
        List<String> list = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            list.add(line);
        }
        fu.close();
        return list;
    }
    
    /**
     * This method is used to pick-up differential expressed genes
     * @throws IOException
     */
    @Test
    public void pickUpDiffExpGenes() throws IOException {
        // Cutoff for p-value
        double pvalueCutOff = 0.01; // To get around 300 genes (5% of total, 6208 genes)
        // File out Control vs Up 
        String controlVsUpFileName = RESULT_DIR + "ControlVsUpTestStat.txt";
        String controlVsDownFileName = RESULT_DIR + "ControlVsDownTestStat.txt";
        // List of genes in up
        List<String> upExpUp = new ArrayList<String>(); // Genes from Up with higher exp
        List<String> upExpLow = new ArrayList<String>(); // Genes from Up with lower exp
        pickUpDiffExpGenes(pvalueCutOff, 
                           controlVsUpFileName, 
                           upExpUp, 
                           upExpLow);
        List<String> downExpLow = new ArrayList<String>(); // Genes from Down with lower exp
        List<String> downExpUp = new ArrayList<String>(); // Genes from Down with higher exp
        pickUpDiffExpGenes(pvalueCutOff, 
                           controlVsDownFileName,
                           downExpUp, 
                           downExpLow);
        // Print out number information
        System.out.println("Diff genes in up samples:");
        System.out.printf("Up: %s; Low: %s%n",
                          upExpUp.size(),
                          upExpLow.size());
        System.out.println("Diff genes in down samples:");
        System.out.printf("Up: %s; low: %s%n",
                           downExpUp.size(),
                           downExpLow.size());
        // Check if any cross-overlapping genes
        // Exp up in up vs exp down in down
        upExpUp.retainAll(downExpLow);
        System.out.println("Exp up in up vs Exp down in down: " + upExpUp.size());
        // Print out these genes
        for (String gene : upExpUp)
            System.out.println("    " + removeQuotes(gene));
        // Exp down in up vs exp up in down
        upExpLow.retainAll(downExpUp);
        System.out.println("Exp down in up vs exp up in down: " + upExpLow.size());
        for (String gene : upExpLow)
            System.out.println("    " + removeQuotes(gene));
        // Make sure these genes are really correct based on UpVsDown test stat
        String upVsDownFileName = RESULT_DIR + "UpVsDownTestStat.txt";
        fu.setInput(upVsDownFileName);
        String line = fu.readLine();
        System.out.println(line);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            if (upExpUp.contains(gene) ||
                upExpLow.contains(gene))
                System.out.println(line);
        }
        fu.close();
    }
    
    private String removeQuotes(String name) {
        return name.substring(1, name.length() - 1);
    }

    private void pickUpDiffExpGenes(double pvalueCutOff,
                                    String controlVsUpFileName,
                                    List<String> upExpUp, List<String> upExpLow)
            throws IOException {
        fu.setInput(controlVsUpFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Use the BY correction.
            double pvalue = Double.parseDouble(tokens[9]);
            if (pvalue > pvalueCutOff)
                continue;
            double tvalue = Double.parseDouble(tokens[1]);
            if (tvalue < 0) // for lower expression
                upExpLow.add(tokens[0]);
            else if (tvalue > 0)
                upExpUp.add(tokens[0]);
        }
        fu.close();
    }
    
    /**
     * This method is used to clean up the filtered samle file so that an output
     * can be loaded into R directly.
     * @throws IOException
     */
    @Test
    public void cleanUpSampleFile() throws IOException {
        String dataFile = RESULT_DIR + "FilteredSampleGeneNormCubicProfile.txt";
        fu.setInput(dataFile);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        List<Integer> neededIndex = new ArrayList<Integer>();
        neededIndex.add(0); // For genes
        int c = 0;
        for (int i = 0; i < headers.length; i++) {
            String header = headers[i];
            if (header.endsWith("AVG_Signal")) {
                neededIndex.add(i);
                System.out.println(c + ": " + header);
                c++;
            }
        }
        FileUtility outFu = new FileUtility();
        String outFileName = RESULT_DIR + "FilteredSampleGeneNormCubicProfile_Avg.txt";
        outFu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        // Handle headers
        for (int i = 0; i < neededIndex.size(); i++) {
            int index = neededIndex.get(i);
            builder.append(headers[index]);
            if (i < neededIndex.size() - 1)
                builder.append("\t");
        }
        outFu.printLine(builder.toString());
        builder.setLength(0);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            for (int i = 0; i < neededIndex.size(); i++) {
                int index = neededIndex.get(i);
                builder.append(tokens[index]);
                if (i < neededIndex.size() - 1)
                    builder.append("\t");
            }
            outFu.printLine(builder.toString());
            builder.setLength(0);
        }
        outFu.close();
        fu.close();
    }
    
    /**
     * This method is used to handle RT-PCR results (created on Dec 13, 2010).
     * @throws Exception
     */
    @Test
    public void generateCountsForValidationDataSet() throws Exception {
        // These samples are used as control: 2FD3, 2F-GFP, 2FD3-siControl and 2F-GFP-siControl
        String[] controls = new String[] {
                "2FD3",
                "2FD3_siCONT",
                "2F-GFP",
                "2F-GFP_siCONT"
        };
        // This sample is used as positive regulated sample: 2F-BRG
        String[] brg1hiSamples = new String[] {
                "2F-BRG1"
        };
        // These samples are used as negative regulated samples: 2F-D3-siBRG1, 2F-GFP-siBRG1, and U3A
        String[] brg1loSamples = new String[] {
                "2FD3-siBRG1",
                "2F-GFP_siBRG1",
        };
        String[] stat1loSample = new String[] {
                "U3A"
        };
        List<String> samples = new ArrayList<String>();
        for (String sample : brg1hiSamples)
            samples.add(sample);
        for (String sample : brg1loSamples)
            samples.add(sample);
        for (String sample : stat1loSample)
            samples.add(sample);
        String dirName = DATA_SET + "Dax/";
        String srcFileName = dirName + "Guanming_RT-qPCR_validationmatrix.txt";
        fu.setInput(srcFileName);
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        List<String> genes = new ArrayList<String>();
        for (int i = 1; i < tokens.length; i++)
            genes.add(tokens[i]);
        System.out.println("Total genes: " + genes.size());
        //Collections.sort(genes);
        //for (String gene : genes)
        //    System.out.println(gene);
        List<Map<String, List<Double>>> geneToSampleToValues = new ArrayList<Map<String, List<Double>>>();
        for (int i = 1; i < tokens.length; i++) 
            geneToSampleToValues.add(new HashMap<String, List<Double>>());
        while ((line = fu.readLine()) != null) {
            tokens =  line.split("\t");
            String sample = tokens[0];
            for (int i = 1; i < tokens.length; i++) {
                Map<String, List<Double>> sampleToValues = geneToSampleToValues.get(i - 1);
                List<Double> values = sampleToValues.get(sample);
                if (values == null) {
                    values = new ArrayList<Double>();
                    sampleToValues.put(sample, values);
                }
                values.add(new Double(tokens[i]));
            }
        }
        fu.close();
        System.out.println("Total samples in values: " + geneToSampleToValues.size());
        String outFileName = dirName + "Guanming_RT-qPCR_validationmatrix_results_with_p_value_cutoff_005.txt";
//        String outFileName = dirName + "Guanming_RT-qPCR_validationmatrix_results_cutoff_005.txt";
        boolean needPvalue = true;
        fu.setOutput(outFileName);
        // Generate a title
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        for (String sample : samples) {
            for (String control : controls) {
                builder.append("\t").append(sample).append(" vs ").append(control);
            }
        }
        fu.printLine(builder.toString());
        double pvalueCutoff = 0.05d;
//        double pvalueCutoff = 0.10d;
        for (int i = 0; i < genes.size(); i++) {
            String gene = genes.get(i);
            builder.setLength(0);
            Map<String, List<Double>> sampleToValues = geneToSampleToValues.get(i);
            builder.append(gene);
            for (String sample : samples) {
                List<Double> sampleValues = sampleToValues.get(sample);
                double[] sampleArray = convertListToArray(sampleValues);
                for (String control : controls) {
                    List<Double> controlValues = sampleToValues.get(control);
                    double[] controlArray = convertListToArray(controlValues);
                    double t = TestUtils.t(sampleArray, controlArray);
                    double pvalue = TestUtils.tTest(sampleArray, controlArray);
                    if (pvalue <= pvalueCutoff) {
                        if (t > 0)
                            builder.append("\t1");
                        else if (t == 0)
                            builder.append("\t0");
                        else 
                            builder.append("\t-1");
                    }
                    else
                        builder.append("\t0");
                    if (needPvalue)
                        builder.append(" (").append(pvalue).append(")");
                }
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    private double[] convertListToArray(List<Double> list) {
        double[] rtn = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            rtn[i] = list.get(i);
        return rtn;
    }
    
    
    /**
     * This method is used to generate files for BioConductor. Three files are needed for
     * BioConductor analysis: two data matrixes, and one gene names. The labels will be created in R as a
     * vector of 0 and 1.
     * @throws IOException
     */
    @Test
    public void generateFilesForBioConductor() throws IOException {
        // These samples are used as control: 2FD3, 2F-GFP, 2FD3-siControl and 2F-GFP-siControl
        String[] controls = new String[] {
                "2FD3",
                "2F-GFP",
                "2FD3_siCONT",
                "2F-GFP_siCONT"
        };
        // This sample is used as positive regulated sample: 2F-BRG
        String[] upSamples = new String[] {
                "2F-BRG1"
        };
        // These samples are used as negative regulated samples: 2F-D3-siBRG1, 2F-GFP-siBRG1, and U3A
        String[] downSamples = new String[] {
                "2FD3-siBRG1",
                "2F-GFP_siBRG1",
                "U3A"
        };
        // Two matrix files should be created, and should be aligned as shown
        String dataFile = RESULT_DIR + "FilteredSampleGeneNormCubicProfile.txt";
        fu.setInput(dataFile);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        Map<String, Map<String, Double>> geneToData = new HashMap<String, Map<String, Double>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // The first one should be gene
            String gene = tokens[0];
            if (gene.length() == 0)
                throw new IllegalStateException("Empty Gene in line: " + line);
            Map<String, Double> values = new HashMap<String, Double>();
            // Check samples
            for (int i = 1; i < tokens.length; i++) {
                // Get header
                String header = headers[i];
                if (header.endsWith(".AVG_Signal")) {
                    values.put(header,
                               Double.parseDouble(tokens[i]));
                }
            }
            geneToData.put(gene, values);
        }
        fu.close();
        // For output
        List<String> geneList = new ArrayList<String>(geneToData.keySet());
        // Do a sorting
        Collections.sort(geneList);
        System.out.println("Total genes: " + geneList.size());
        // Gene lists
        String geneNameList = RESULT_DIR + "GeneList.txt";
        fu.setOutput(geneNameList);
        for (String gene : geneList)
            fu.printLine(gene);
        fu.close();
        // Get the first matrix: control vs up-regulated
        String profileName = RESULT_DIR + "ControlVsUpProfile.txt";
        String labelFileName = RESULT_DIR + "ControlVsUpLabels.txt";
        generateControlVsSampleFiles(controls, 
                                     upSamples,
                                     headers,
                                     geneToData,
                                     geneList,
                                     profileName,
                                     labelFileName);
        // Get the second matrix: control vs down-regulated
        profileName = RESULT_DIR + "ControlVsDownProfile.txt";
        labelFileName = RESULT_DIR + "ControlVsDownLabels.txt";
        generateControlVsSampleFiles(controls, 
                                     downSamples, 
                                     headers, 
                                     geneToData, 
                                     geneList, 
                                     profileName, 
                                     labelFileName);
        profileName = RESULT_DIR + "UpVsDownProfile.txt";
        labelFileName = RESULT_DIR + "UpVsDownLabels.txt";
        generateControlVsSampleFiles(upSamples, 
                                     downSamples, 
                                     headers, 
                                     geneToData, 
                                     geneList, 
                                     profileName, 
                                     labelFileName);        
    }

    private void generateControlVsSampleFiles(String[] controls,
                                              String[] samples,
                                              String[] headers,
                                              Map<String, Map<String, Double>> geneToData,
                                              List<String> geneList,
                                              String outFileName,
                                              String labelFileName) throws IOException {
        fu.setOutput(outFileName);
        List<Double> valueList = new ArrayList<Double>();
        StringBuilder builder = new StringBuilder();
        for (String gene : geneList) {
            Map<String, Double> headerToValue = geneToData.get(gene);
            valueList.clear();
            // Control values
            for (String control : controls) {
                String headerLabel = control + ".AVG_Signal";
                // Search through values
                // Use the headers array to ensure the correct order
                for (String header : headers) {
                    Double value = headerToValue.get(header);
                    if (value == null)
                        continue;
                    if (header.endsWith(headerLabel))
                        valueList.add(value);
                }
            }
            // up samples
            for (String up : samples) {
                String headerLabel = up + ".AVG_Signal";
                for (String header : headers) {
                    Double value = headerToValue.get(header);
                    if (value == null)
                        continue;
                    if (header.endsWith(headerLabel))
                        valueList.add(value);
                }
            }
            // Output value list
            builder.setLength(0);
            for (Iterator<Double> it = valueList.iterator(); it.hasNext();) {
                Double value = it.next();
                builder.append(value);
                if (it.hasNext())
                    builder.append("\t");
            }
            fu.printLine(builder.toString());
        }
        fu.close();
        // Output labels
        fu.setOutput(labelFileName);
        builder.setLength(0);
        // All samples have been tripled
        for (String control : controls) {
            for (int i = 0; i < 3; i++)
                builder.append(0).append("\t");
        }
        for (int i = 0; i < samples.length; i++) {
            if (i < samples.length - 1) {
                for (int j = 0; j < 3; j++) {
                    builder.append(1).append("\t");
                }
            }
            else {
                builder.append(1).append("\t").append(1).append("\t").append(1);
            }
        }
        fu.printLine(builder.toString());
        fu.close();
    }
    
    /**
     * This method is used to filter the data set based on detection pvalue. Only genes
     * with p-value < 0.01 are used in this file.
     * @throws IOException
     */
    @Test
    public void filterDatasets() throws IOException {
        String outFileName = RESULT_DIR + "FilteredSampleGeneNormCubicProfile.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        fu.setInput(SOURCE_FILE);
        String line = fu.readLine();
        outFu.printLine(line);
        String[] headers = line.split("\t");
        boolean isGood = false;
        int totalRow = 0;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].length() == 0)
                continue; // No gene mapped
            isGood = true;
            // Check p-value
            for (int i = 0; i < tokens.length; i++) {
                if (headers[i].endsWith("Detection Pval")) {
                    // P-value
                    double pvalue = Double.parseDouble(tokens[i]);
                    if (pvalue > 0.01) {
                        isGood = false;
                        break;
                    }
                }
            }
            if (isGood) {
                outFu.printLine(line);
                totalRow ++;
                isGood = false;
            }
        }
        outFu.close();
        fu.close();
    }
    
    @Test
    public void checkNumberOfUsableGenes() throws IOException {
        //fu.setInput(SOURCE_FILE);
        //String fileName = DATA_SET  + "mohamed/Avg_Norm_Bg_Sub.txt";
        String fileName = DATA_SET + "Dax/avg_norm_bg_sub.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get the headers
        String[] headers = line.split("\t");
        Set<String> genes = new HashSet<String>();
        boolean isGood = false;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].length() == 0)
                continue; // No gene mapped
            isGood = true;
            // Check p-value
            for (int i = 0; i < tokens.length; i++) {
                if (headers[i].endsWith("Detection Pval")) {
                    // P-value
                    double pvalue = Double.parseDouble(tokens[i]);
                    if (pvalue > 0.05) {
                        isGood = false;
                        break;
                    }
                }
            }
            if (isGood) {
                genes.add(tokens[0]);
                isGood = false;
            }
        }
        fu.close();
        System.out.print("Total usable genes: " + genes.size());
    }
    
    @Test
    public void checkHeaders() throws IOException {
        fu.setInput(SOURCE_FILE);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        for (String header : headers)
            System.out.println(header);
        fu.close();
    }
    
    private class ValidationResult implements StatisticalSummary {

        double mean;
        double sd;
        
        
        public double getMax() {
            
            return 0;
        }

        public double getMean() {
            return mean;
        }

        public double getMin() {
            
            return 0;
        }

        public long getN() {
            return 3;
        }

        public double getStandardDeviation() {
            return sd;
        }

        public double getSum() {
            
            return 0;
        }

        public double getVariance() {
            return sd * sd;
        }
    }
}
