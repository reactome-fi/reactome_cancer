/*
 * Created on Nov 18, 2009
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.util.FileUtility;

/**
 * This class is used to process jobs related to R.
 * @author wgm
 *
 */
public class RUtilities {
    private static FileUtility fu = new FileUtility();
    
    public static Map<String, Double> loadGeneToValue(String fileName,
                                                      int valueIndex) throws IOException {
        fu.setInput(fileName);
        String line = null;
        Map<String, Double> rtn = new HashMap<String, Double>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            rtn.put(tokens[0],
                    new Double(tokens[valueIndex]));
        }
        return rtn;
    }
    
    public static List<String> loadGeneList(String fileName) throws IOException {
        fu.setInput(fileName);
        String line = fu.readLine();
        String[] genes = line.split(", ");
        fu.close();
        List<String> geneList = Arrays.asList(genes);
        return new ArrayList<String>(geneList); // Use array list for common operations
    }
    
    public static List<Set<String>> loadClustersFromDynamicTOM(String clusterFile,
                                                               String matrixFile) throws IOException {
        // Load all genes
        fu.setInput(matrixFile);
        String line = fu.readLine();
        fu.close();
        // Get gene list
        String[] tokens = line.split("\t");
        List<String> geneList = new ArrayList<String>();
        for (String token : tokens) {
            if (token.length() == 0)
                continue;
            geneList.add(token);
        }
        System.out.println("Total genes: " + geneList.size());
        // Get the clusters
        fu.setInput(clusterFile);
        line = fu.readLine();
        fu.close();
        tokens = line.split(", ");
        List<Integer> clusters = new ArrayList<Integer>();
        int totalCluster = 0;
        for (String token : tokens) {
            int cluster = new Integer(token);
            if (cluster > totalCluster)
                totalCluster = cluster;
            clusters.add(cluster);
        }
        System.out.println("Total cluster numbers: " + clusters.size());
        List<Set<String>> rtn = new ArrayList<Set<String>>();
        for (int i = 0; i < totalCluster; i++)
            rtn.add(new HashSet<String>());
        for (int i = 0; i < geneList.size(); i++) {
            String gene = geneList.get(i);
            Integer cluster = clusters.get(i);
            if (cluster == 0)
                continue; // Not in any cluster
            // Remove -1
            Set<String> geneCluster = rtn.get(cluster - 1);
            geneCluster.add(gene);
        }
        return rtn;
    }
    
    /**
     * Lost a list of samples output from R hclust.
     * @param fileName
     * @return
     * @throws IOException
     */
    public static List<String> loadSamplesInRHClust(String fileName) throws IOException {
        fu.setInput(fileName);
        List<String> samples = new ArrayList<String>();
        String line = null;
        while ((line = fu.readLine()) != null) {
            line = line.trim();
            samples.add(line);
        }
        fu.close();
        return samples;
    }
    
    /**
     * This method is used to load cluster files dumped out from R hclust.
     * @param fileName
     * @return
     * @throws IOException
     */
    public static List<List<String>> loadHierarchicalClustersFromR(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        List<List<String>> clusters = new ArrayList<List<String>>();
        fu.setInput(fileName);
        String line = null;
        List<String> cluster = null;
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                continue;
            String[] tokens = line.split("\t");
            cluster = new ArrayList<String>();
            for (String token : tokens)
                cluster.add(token);
            clusters.add(cluster);
        }
        fu.close();
        return clusters;
    }
}
