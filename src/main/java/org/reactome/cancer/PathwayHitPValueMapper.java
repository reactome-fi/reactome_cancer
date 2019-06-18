/*
 * Created on Jun 8, 2007
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.reactome.r3.util.FileUtility;

/**
 * This class is used to map the calculated p-values from binomial test to p-values
 * generated from random test.
 * @author guanming
 *
 */
public class PathwayHitPValueMapper {
    private double cutoff;
    private double minimum;
    // Tested threshold values used as steps
    private List<Double> thresholdValues;
    private Map<String, List<MappedValue>> topicToMappedValues;
    
    public PathwayHitPValueMapper() {
        topicToMappedValues = new HashMap<String, List<MappedValue>>();
        thresholdValues = new ArrayList<Double>();
    }
    
    public void setMinimum(double min) {
        this.minimum = min;
    }
    
    public void setCutOffValue(double cutoff) {
        this.cutoff = cutoff;
    }
    
    /**
     * The returned type is String so that less than minimum can be returned.
     * @param pathway
     * @param pValue
     * @return
     */
    public String map(String pathway,
                      double pValue) {
        // If the passed pValue is greater than cutoff,
        // just return 1.0. Not interesting!
        if (pValue > cutoff)
            return 1.0 + "";
        List<MappedValue> mappedValues = topicToMappedValues.get(pathway);
        if (mappedValues == null || mappedValues.size() == 0)
            return "<" + minimum;
        // Find the closed pValue from the threshold list
        double diff = Double.MAX_VALUE;
        double closedThreshold = -1.0d;
        for (Double d : thresholdValues) {
            double tmpDiff = d - pValue;
            if (tmpDiff > 0 && tmpDiff < diff) {
                diff = tmpDiff;
                closedThreshold = d;
            }
        }
        // Need to find the closest threshold that is less than the passed pValue
        MappedValue closedValue = null;
        // pValue should be less than threshold. Diff should be positive always
        for (int i = 0; i < mappedValues.size(); i++) {
            MappedValue tmp = mappedValues.get(i);
            if (tmp.threshold == closedThreshold) {
                closedValue = tmp;
                break;
            }
        }
        if (closedValue == null) {
            // meaning no occurences at the closed threshold
            return "<" + minimum;
        }
        if (diff > 0 && closedValue != null)
            return closedValue.pValue + "";
        // Should not be here
        return null;
    }
    
    public void loadData(String fileName) throws IOException {
        topicToMappedValues.clear();
        thresholdValues.clear();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int permutationNumber = 0;
        double threshold = 0.0d;
        // used to extract out threshold
        int index1, index2;
        // This RegExp is used to get the permutation number
        String regExp = "(\\d+)$";
        Pattern pattern = Pattern.compile(regExp);
        while ((line = fu.readLine()) != null) {
            if (line.equals("------"))
                break;
            if (line.length() == 0 ||
                line.startsWith("Total Topics"))
                continue;
            if (line.startsWith("p-value:")) {
                // Need to get permutation number
                Matcher matcher = pattern.matcher(line);
                // Have to call find() method first
                matcher.find();
                String permutation = matcher.group(1);
                permutationNumber = Integer.parseInt(permutation);
                index1 = line.indexOf(":");
                index2 = line.indexOf(",");
                threshold = Double.parseDouble(line.substring(index1 + 1, index2).trim());
                thresholdValues.add(threshold);
            }
            else {
                String[] tokens = line.split("\t");
                // Translation(R)   7   119 0.018772677078403535
                double ratio = Double.parseDouble(tokens[1]) / permutationNumber;
                List<MappedValue> values = topicToMappedValues.get(tokens[0]);
                if (values == null) {
                    values = new ArrayList<MappedValue>();
                    topicToMappedValues.put(tokens[0], values);
                }
                MappedValue mappedValue = new MappedValue();
                mappedValue.threshold = threshold;
                mappedValue.pValue = ratio;
                values.add(mappedValue);
            }
        }
        fu.close();
    }
    
    private class MappedValue {
        private double threshold;
        private double pValue;
    }
    
}
