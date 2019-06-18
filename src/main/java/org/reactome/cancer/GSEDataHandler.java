/*
 * Created on Feb 7, 2011
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.GeneExpressionDataSet;

/**
 * A generalization for handling the GSE data sets.
 * @author wgm
 *
 */
public abstract class GSEDataHandler {
    
    protected final FileUtility fu = new FileUtility();

    public GSEDataHandler() {
        
    }
    
    /**
     * This method is used to map probeset to genes for a passed GeneExpressionDataSet object.
     * @param dataset
     * @param probesetToGene
     * @throws IOException
     */
    public void mapProbesetToGenes(GeneExpressionDataSet dataset,
                                   Map<String, String> probesetToGene) throws IOException {
        System.out.println("Features in data set: " + dataset.getFeatureList().size());
        System.out.println("Values in data set: " + dataset.getValues().size());
        List<String> featureList = dataset.getFeatureList();
        // Need to copy to avoid ConcurrentModificationException
        List<String> copy = new ArrayList<String>(featureList);
        for (String probeset : copy) {
            String gene = probesetToGene.get(probeset);
            if (gene == null)
                dataset.dropFeature(probeset);
            else {
                dataset.replaceFeatureName(probeset, gene);
            }
        }
    }
    
    /**
     * After probesets are mapped to genes, some genes may appear multiple times in the dataset.
     * In this method, expression values for these genes will be averaged.
     * @param dataset
     */
    public GeneExpressionDataSet averageValuesForSameGenes(GeneExpressionDataSet dataset) {
        List<String> featureList = dataset.getFeatureList();
        System.out.println("Total featureList before averaging: " + featureList.size());
        Map<String, List<Integer>> featureToIndex = new HashMap<String, List<Integer>>();
        for (int i = 0; i < featureList.size(); i++) {
            String feature = featureList.get(i);
            List<Integer> list = featureToIndex.get(feature);
            if (list == null) {
                list = new ArrayList<Integer>();
                featureToIndex.put(feature, list);
            }
            list.add(i);
        }
        // A non-redundant gene list
        List<String> geneList = new ArrayList<String>(featureToIndex.keySet());
        Collections.sort(geneList);
        System.out.println("Total genes: " + geneList.size());
        GeneExpressionDataSet rtn = new GeneExpressionDataSet();
        rtn.setSampleList(dataset.getSampleList());
        rtn.setFeatureList(geneList);
        for (String gene : geneList) {
            if (gene.equals("GDA"))
                System.out.println();
            List<Integer> indices = featureToIndex.get(gene);
            List<Double> average = averageValues(indices, dataset);
            rtn.addValuesForFeature(average);
        }
        return rtn;
    }
    
    private List<Double> averageValues(List<Integer> indices,
                                       GeneExpressionDataSet dataset) {
        List<Double> rtn = new ArrayList<Double>();
        int colSize = dataset.getSampleList().size();
        double total = 0.0;
        int count = 0;
        for (int i = 0; i < colSize; i++) {
            total = 0.0;
            count = 0;
            for (int j = 0; j < indices.size(); j++) {
                List<Double> original = dataset.getValuesForFeature(indices.get(j));
//                if (i == original.size())
//                    System.out.println("Wrong!");
                Double value = original.get(i);
                if (value != null) {
                    total += value;
                    count ++;
                }
            }
            if (count > 0)
                rtn.add(total / count);
            else
                rtn.add(null);
        }
        return rtn;
    }
    
    public abstract Map<String, String> loadProbeIdToGene(String fileName) throws IOException; 
    
}
