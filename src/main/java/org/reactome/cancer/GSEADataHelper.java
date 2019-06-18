/*
 * Created on Feb 15, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to process GSEA related data and files.
 * @author wgm
 *
 */
public class GSEADataHelper {
    private FileUtility fu = new FileUtility();
    
    /**
     * This method is used to generate a class file for GSEA.
     */
    @Test
    public void generateClassFile() throws IOException {
        //String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_110609/";
        //String expFileName = dirName + "TCGA_17-19_median_excluded_controls.txt";
        String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_072610/";
        String expFileName = dirName + "TCGA_batch9-15_17-19_21-22_24.UE.txt";
        fu.setInput(expFileName);
        String line = fu.readLine();
        fu.close();
        String[] tokens = line.split("\t");
        System.out.println("Total samples: " + (tokens.length - 1));
        //String outFileName = dirName + "TCGA_17-19_median_excluded_controls.cls";
        String outFileName = dirName + "TCGA_batch9-15_17-19_21-22_24.UE.cls";
        fu.setOutput(outFileName);
        fu.printLine((tokens.length - 1) + " 2 1");
        fu.printLine("# NORMAL OV");
        StringBuilder builder = new StringBuilder();
        for (int i = 1; i < tokens.length; i++) {
            String name = tokens[i];
            //if (name.endsWith("-11")) {
            if (name.contains(".11A.")) {
                builder.append("0 ");
            }
            else
                builder.append("1 ");
        }
        fu.printLine(builder.toString());
        fu.close();
    }
    
    /**
     * This method is used to remove NA rows. GSEA cannot handle NA.
     * @throws IOException
     */
    @Test
    public void removeNAForGSEA() throws IOException {
        String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_110609/";
//        String expFileName = dirName + "TCGA_Batch9-15_median_Symbol_excluded.txt";
//        String outFileName = dirName + "TCGA_Batch9-15_median_Symbol_excluded_na_removed.txt";
        String expFileName = dirName + "TCGA_17-19_median_excluded.txt";
        String outFileName = dirName + "TCGA_17-19_median_excluded_na_removed.txt";
        fu.setInput(expFileName);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        String line = fu.readLine();
        outFu.printLine(line);
        int index = 0;
        String subLine;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            subLine = line.substring(index + 1);
            if (!subLine.contains("NA"))
                outFu.printLine(line);
        }
        outFu.close();
        fu.close();
    }
    
    /**
     * This method is used to copy gene expression for controls to another files in order to create
     * a new gene expression for use by GSEA.
     * @throws IOException
     */
    @Test
    public void mergeControlsToSamplesForGeneExp() throws IOException {
        String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_110609/";
        String expFileName = dirName + "TCGA_Batch9-15_median_Symbol_excluded.txt";
        // Get a list of samples
        fu.setInput(expFileName);
        String line = fu.readLine();
        fu.close();
        // Create a list of samples for output
        List<String> samples = new ArrayList<String>();
        String[] tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++) {
            if (tokens[i].endsWith("-11"))
                samples.add(tokens[i]);
        }
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> oldGeneToExp = helper.loadGeneExp(expFileName);
        expFileName = dirName + "TCGA_17-19_median_excluded.txt";
        fu.setInput(expFileName);
        line = fu.readLine();
        fu.close();
        tokens = line.split("\t");
        for (int i = 1; i < tokens.length; i++) {
            samples.add(tokens[i]);
        }
        Map<String, Map<String, Double>> newGeneToExp = helper.loadGeneExp(expFileName);
        List<String> genes = new ArrayList<String>(newGeneToExp.keySet());
        genes.retainAll(oldGeneToExp.keySet());
        Collections.sort(genes);
        
        String outFileName = dirName + "TCGA_17-19_median_excluded_controls.txt";
        fu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        for (String sample : samples)
            builder.append("\t").append(sample);
        fu.printLine(builder.toString());
        builder.setLength(0);
        for (String gene : genes) {
            Map<String, Double> oldExp = oldGeneToExp.get(gene);
            Map<String, Double> newExp = newGeneToExp.get(gene);
            builder.append(gene);
            for (String sample : samples) {
                Double value = oldExp.get(sample);
                if (value == null)
                    value = newExp.get(sample);
                builder.append("\t").append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * Use this method to do row-based median centered transformation for a set of gene expression.
     * Note: the following method is not correct. The method to do median centering is wrong.
     * @throws IOException
     */
    @Test
    public void globalBasedMedianCenteredGeneExp() throws IOException {
//        String dirName = "datasets/TCGA/OvarianCancer/TothillDataset/";
//        String srcFileName = dirName + "Tothill_TCGA_plus2.txt";
//        fu.setInput(srcFileName);
//        FileUtility outFu = new FileUtility();
//        String outFileName = dirName + "Tothill_TCGA_plus2_global_median_centered.txt";
//        outFu.setOutput(outFileName);
//        String line = fu.readLine();
//        outFu.printLine(line);
//        DescriptiveStatistics stat = new DescriptiveStatistics();
//        StringBuilder builder = new StringBuilder();
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            for (int i = 1; i < tokens.length; i++) {
//                stat.addValue(Math.log(new Double(tokens[i])));
//            }
//        }
//        fu.close();
//        double median = stat.getPercentile(0.50d);
//        fu.setInput(srcFileName);
//        line = fu.readLine();
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            builder.setLength(0);
//            builder.append(tokens[0]);
//            for (int i = 1; i < tokens.length; i++) {
//                double value = Double.parseDouble(tokens[i]);
//                value = (value - median) / median;
//                builder.append("\t").append(value);
//            }
//            outFu.printLine(builder.toString());
//        }
//        fu.close();
//        outFu.close();
    }
    
}
