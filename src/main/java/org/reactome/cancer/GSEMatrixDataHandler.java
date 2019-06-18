/*
 * Created on Jan 21, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.reactome.r3.util.GeneExpressionDataSet;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to handle GSE matrix data files downloaded from the GEO data repository.
 * @author wgm
 *
 */
public class GSEMatrixDataHandler extends GSEDataHandler {
    public GSEMatrixDataHandler() {
    }
    
    /**
     * Use this method to actually run the load, map, and average.
     * @throws IOException
     */
    @Test
    public void processMatrixFiles() throws IOException {
        String dirName = "/Users/wgm/datasets/BreastCancer/GSE3143/";
        String gseName = dirName + "GSE3143_series_matrix.txt";
        int headerRow = 50;
        int finalRow = 12630;
        String gplName = dirName + "GPL8300-39822.txt";
//        String outName = dirName + "GSE3143_MappedGenes_z.txt";
//        String outName = dirName + "GSE3143_MappedGenes_log_z.txt";
        String outName = dirName + "GSE3143_MappedGenes_z_on_samples_012111.txt";
//        String outName = dirName + "GSE3143_MappedGenes_log_z_on_samples_012111.txt";
        
        GeneExpressionDataSet dataset = loadGeneExpressionDataSet(gseName, 
                                                                  finalRow,
                                                                  headerRow);
        Map<String, String> probesetToGene = loadProbeIdToGene(gplName);
        mapProbesetToGenes(dataset, probesetToGene);
        GeneExpressionDataSet averagedDataset = averageValuesForSameGenes(dataset);
//        averagedDataset.logTransformation();
//        averagedDataset.zscoreTansformation();
        averagedDataset.sampleWiseZscoreTransformation();
        averagedDataset.export(outName);
    }
    
    @Override
    public Map<String, String> loadProbeIdToGene(String fileName) throws IOException {
        fu.setInput(fileName);
        Map<String, String> probesetToGenes = new HashMap<String, String>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue; // Escape comment and annotation lines
            String[] tokens = line.split("\t");
            if (tokens.length < 10)
                continue;
            String gene = tokens[10].trim();
            if (gene.length() == 0 || gene.contains("///")) 
                continue; // Don't need to pick up unmapped probesets and multiple mapped genes
            probesetToGenes.put(tokens[0], gene);
        }
        fu.close();
        return probesetToGenes;
    }

    /**
     * Load a geo matrix file as a GeneExpressionDataset object.
     * @param fileName the matrix file to be loaded
     * @param finalRow the last row for loading. This row should NOT be loaded (exclusively)
     * @param headRow the row contains sample information as header.
     * @return
     * @throws IOException
     */
    public GeneExpressionDataSet loadGeneExpressionDataSet(String fileName,
                                                           int finalRow,
                                                           int headRow) throws IOException {
        fu.setInput(fileName);
        GeneExpressionDataSet dataset = new GeneExpressionDataSet();
        String line = null;
        int rowCount = 0;
        while ((line = fu.readLine()) != null) {
            rowCount ++;
            line = line.trim();
            if (rowCount == finalRow)
                break;
            if (line.startsWith("!") || line.length() == 0) {
                continue;
            }
            String[] tokens = line.split("\t");
            if (rowCount == headRow) {
                // Get samples
                for (int i = 1; i < tokens.length; i++) {
                    String sample = tokens[i];
                    sample = InteractionUtilities.removeQuotationMarks(sample);
                    dataset.addSample(sample);
                }
            }
            else {
                // Actual values
                String feature = InteractionUtilities.removeQuotationMarks(tokens[0]);
                dataset.addFeature(feature);
                List<Double> featureValues = new ArrayList<Double>();
                for (int i = 1; i < tokens.length; i++) {
                    String value = tokens[i];
                    if (value.length() == 0 || value.equals("null"))
                        featureValues.add(null);
//                    else if (value.equals("0")) {
//                        System.out.println("Value is zero: " + feature + ", " + feature);
//                        featureValues.add(null);
//                    }
                    else
                        featureValues.add(new Double(value));
                }
                dataset.addValuesForFeature(featureValues);
            }
        }
        fu.close();
        return dataset;
    }
    
    public GeneExpressionDataSet loadGeneExpressionDataSet(String fileName) throws IOException {
        // Find the first and the last row
        int firstRow = 0;
        int lastRow = 0;
        fu.setInput(fileName);
        String line = null;
        int rowIndex = 0;
        while ((line = fu.readLine()) != null) {
            rowIndex ++;
            if (line.startsWith("\"ID_REF\""))
                firstRow = rowIndex;
            else if (line.startsWith("!series_matrix_table_end"))
                lastRow = rowIndex;
        }
        fu.close();
        return loadGeneExpressionDataSet(fileName, lastRow, firstRow);
    }
    
}
