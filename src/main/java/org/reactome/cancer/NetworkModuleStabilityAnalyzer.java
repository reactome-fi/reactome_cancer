/*
 * Created on Jan 24, 2013
 *
 */
package org.reactome.cancer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.junit.Test;
import org.reactome.annotate.AnnotationHelper;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.r3.graph.SpectralPartitionNetworkCluster;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to test some module stability as requested by Lincoln. There are
 * two types of issues here: pathway annotations and module stability. For pathway annotation
 * analysis, pathways have FDR less than 0.001 are selected from any selected modules, and then
 * a Hypogeometric test is done. For modules, module pair-wise comparison between two lists
 * are done. In order to select modules for analysis, module size, an adjusted parameter, is
 * used to select modules.
 * @author gwu
 *
 */
public class NetworkModuleStabilityAnalyzer {
    private final int MODULE_SIZE_FILTER = 5;
    private final int SAMPLE_SIZE_FILTER = 4;
    private final double P_VALUE_CUTOFF = 0.001d;
    private final double ANNOTATION_FDR_VALUE_CUTOFF = 0.0001d;
    private FileUtility fu = new FileUtility();
    
    public NetworkModuleStabilityAnalyzer() {
    }
    
    private Map<String, Set<String>> generateRandomSamples(Map<String, Set<String>> realSampleToGenes,
                                                           Set<String> allGenes) {
        Map<String, Set<String>> randomSampleToGenes = new HashMap<String, Set<String>>();
        RandomDataImpl randomizer = new RandomDataImpl();
        for (String sample : realSampleToGenes.keySet()) {
            Set<String> realGenes = realSampleToGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(realGenes, allGenes);
            if (shared.size() == 0)
                randomSampleToGenes.put(sample, new HashSet<String>());
            else {
                Set<String> randomGenes = MathUtilities.randomSampling(allGenes,
                                                                       shared.size(),
                                                                       randomizer);
                randomSampleToGenes.put(sample, randomGenes);
            }
        }
        return randomSampleToGenes;
    }
    
    @Test
    public void testAnnotation() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes();
        
        List<Set<String>> clusters = cluster(fis, sampleToGenes);
        System.out.println("Number of clusters: " + clusters.size());
        
        Set<String> pathways = annotateClusters(clusters);
        System.out.println("Total pathways: " + pathways.size());
    }
    
    private Set<String> annotateClusters(List<Set<String>> clusters) throws Exception {
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        Set<String> pathways = new HashSet<String>();
        for (Set<String> cluster : clusters) {
            List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(cluster,
                                                                                 AnnotationType.Pathway);
            for (GeneSetAnnotation annotation : annotations) {
                String fdr = annotation.getFdr();
                if (fdr.startsWith("<"))
                    pathways.add(annotation.getTopic());
//                else if (new Double(fdr) < ANNOTATION_FDR_VALUE_CUTOFF)
//                    pathways.add(annotation.getTopic());
            }
        }
        return pathways;
    }
    
    private Set<String> getSelectedGenes(Map<String, Set<String>> sampleToGenes,
                                         Set<String> totalGenes) {
        Set<String> genesInSamples = InteractionUtilities.grepAllGenes(sampleToGenes);
        return InteractionUtilities.getShared(genesInSamples, totalGenes);
    }
    
    @Test
    public void checkModuleStability() throws Exception {
        // The following statement is used to get the total number of pathways
        AnnotationHelper annotationHelper = new AnnotationHelper();
        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToTermsMap(AnnotationType.Pathway);
        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
        int totalPathway = pathwayToGenes.size();
        System.out.println("Total pathways: " + totalPathway);
        
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes();
        
        // Used for a random sample set
        sampleToGenes = generateRandomSamples(sampleToGenes, totalGenes);
        Set<String> selectedGenes = getSelectedGenes(sampleToGenes, totalGenes);
        List<Set<String>> clusters = cluster(fis, sampleToGenes);
        Set<String> pathways = annotateClusters(clusters);
        System.out.println("Number of clusters: " + clusters.size());
        System.out.println("Total pathways: " + pathways.size());
//        if (true)
//            return;
        // Want to do a sample based random sampling
        RandomDataImpl randomizer = new RandomDataImpl();
        double[] percentages = new double[] {0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20};
        int permutation = 10;
        SummaryStatistics moduleStats = new SummaryStatistics();
        SummaryStatistics sharedModuleStats = new SummaryStatistics();
        SummaryStatistics sharedModulePvalueStats = new SummaryStatistics();
        SummaryStatistics pathwayStats = new SummaryStatistics();
        SummaryStatistics pathwaySharedStats = new SummaryStatistics();
        SummaryStatistics pathwayPValueStats = new SummaryStatistics();
        SummaryStatistics geneStats = new SummaryStatistics();
        StringBuilder builder = new StringBuilder();
        builder.append("Percentage\tSample_Size\t" +
                "Selected_Genes\tSD\t" + 
                "Modules\tSD\t" + 
                "Shared_Module\tSD\t" +
                "pValue_Mean\tSD\t" +
                "Annotations\tSD\t" + 
                "Shared_Annotations\tSD\t" +
                "pValue_Mean\tSD\n");
        builder.append("1.00\t").append(sampleToGenes.size()).append("\t");
        builder.append(selectedGenes.size()).append("\t0\t");
        builder.append(clusters.size()).append("\t0\t");
        builder.append("1\t0\t");
        builder.append("0\t0\t");
        builder.append(pathways.size()).append("\t0\t");
        builder.append(pathways.size()).append("\t0\t");
        builder.append("0\t0\n");
        for (double percentage : percentages) {
            System.out.println("\nSample percentage: " + percentage);
            int sampleSize = (int)(sampleToGenes.size() * percentage);
            moduleStats.clear();
            sharedModuleStats.clear();
            sharedModulePvalueStats.clear();
            pathwayStats.clear();
            pathwaySharedStats.clear();
            pathwayPValueStats.clear();
            geneStats.clear();
            for (int i = 0; i < permutation; i++) {
                Set<String> samples1 = MathUtilities.randomSampling(sampleToGenes.keySet(),
                                                                    sampleSize,
                                                                    randomizer);
                Map<String, Set<String>> sampleToGenes1 = new HashMap<String, Set<String>>();
                for (String sample : samples1) {
                    sampleToGenes1.put(sample, sampleToGenes.get(sample));
                }
                // Gene sharing without considering clusters
                Set<String> selectedGenes1 = getSelectedGenes(sampleToGenes1, totalGenes);
                geneStats.addValue(selectedGenes1.size());
                // Clusters sharing
                List<Set<String>> clusters1 = cluster(fis, sampleToGenes1);
                moduleStats.addValue(clusters1.size());
                double[] shared = calculateSimilarModules(clusters,
                                                          clusters1,
                                                          totalGenes.size());
                sharedModuleStats.addValue(shared[0]);
                sharedModulePvalueStats.addValue(shared[1]);
                // Annotation sharing
                Set<String> pathways1 = annotateClusters(clusters1);
                double[] pathwayShared = calculateSimilarAnnotations(pathways,
                                                                     pathways1,
                                                                     totalPathway);
                pathwayStats.addValue(pathways1.size());
                pathwaySharedStats.addValue(pathwayShared[0]);
                pathwayPValueStats.addValue(pathwayShared[1]);
            }
            builder.append(percentage + "\t" + sampleSize + "\t" + 
                    geneStats.getMean() + "\t" + geneStats.getStandardDeviation() + "\t" +
                    moduleStats.getMean() + "\t" + moduleStats.getStandardDeviation() + "\t" +
                    sharedModuleStats.getMean() + "\t" + sharedModuleStats.getStandardDeviation() + "\t" + 
                    sharedModulePvalueStats.getMean() + "\t" + sharedModulePvalueStats.getStandardDeviation() + "\t" + 
                    pathwayStats.getMean() + "\t" + pathwayStats.getStandardDeviation() + "\t" +
                    pathwaySharedStats.getMean() + "\t" + pathwaySharedStats.getStandardDeviation() + "\t" + 
                    pathwayPValueStats.getMean() + "\t" + pathwayPValueStats.getStandardDeviation() + "\n");
        }
        System.out.println(builder.toString());
    }
    
    private double[] calculateSimilarAnnotations(Set<String> refAnnotations,
                                                 Set<String> annotations,
                                                 int totalPathway) throws MathException {
        Set<String> shared = InteractionUtilities.getShared(refAnnotations, annotations);
        double pvalue = MathUtilities.calculateHypergeometricPValue(totalPathway,
                                                                    refAnnotations.size(),
                                                                    annotations.size(),
                                                                    shared.size());
        return new double[]{shared.size(), pvalue};
    }
    
    private double[] calculateSimilarModules(List<Set<String>> refModules,
                                             List<Set<String>> modules,
                                             int totalGene) throws org.apache.commons.math.MathException {
        int count = 0;
        // Used to calculate average
        double logPvalueTotal = 0.0d;
        for (Set<String> refModule : refModules) {
            List<Double> pvalues = new ArrayList<Double>();
            for (Set<String> module : modules) {
                Set<String> shared = InteractionUtilities.getShared(refModule, module);
                double pvalue = MathUtilities.calculateHypergeometricPValue(totalGene,
                                                                            refModule.size(),
                                                                            module.size(),
                                                                            shared.size());
                pvalues.add(pvalue);
            }
            Collections.sort(pvalues);
            Double pvalue = null;
            if (pvalues.size() > 0)
                pvalue = pvalues.get(0);
            else
                pvalue = 1.0d;
            logPvalueTotal += Math.log(pvalue);
            if (pvalue < P_VALUE_CUTOFF) {
                count ++;
            }
        }
        return new double[]{(double) count / refModules.size(),
                Math.exp(logPvalueTotal / refModules.size())};
    }
    
    private List<Set<String>> cluster(Set<String> fis,
                                      Map<String, Set<String>> sampleToGenes) {
        System.out.println("Total samples: " + sampleToGenes.size());
        Set<String> genes = CancerAnalysisUtilitites.selectGenesInSamples(SAMPLE_SIZE_FILTER, 
                                                                          sampleToGenes);
        System.out.println("Total genes: " + genes.size());
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        System.out.println("Total fis: " + fisInGenes.size());
        
        if (fisInGenes.size() == 0)
            return new ArrayList<Set<String>>();
        
        SpectralPartitionNetworkCluster clustering = new SpectralPartitionNetworkCluster();
        List<Set<String>> clusters = clustering.cluster(fisInGenes);
        
        NetworkClusterAnalyzer clusterAnalyzer = new NetworkClusterAnalyzer();
        clusterAnalyzer.filterNetworkClusters(clusters, MODULE_SIZE_FILTER);
        return clusters;
    }
    
    private Map<String, Set<String>> loadSampleToMutatedGenes() throws Exception {
        // For TCGA OV
//        String dirName = R3Constants.OVARIAN_DIR_NAME + "FireHose/";
//        String fileName = dirName + "gdac.broadinstitute.org_OV.Mutation_Significance.Level_4.2012072500.0.0/OV.final_analysis_set.maf";
        // TCGA BRCA
        String dirName = "/Users/gwu/datasets/TCGA/BRCA/gdac.broadinstitute.org_BRCA-TP.Mutation_Assessor.Level_4.2012122100.0.0/";
        String fileName = dirName + "BRCA-TP.maf.annotated";
        
        MATFileLoader fileLoader = new MATFileLoader();
        return fileLoader.loadSampleToGenes(fileName, false);
    }
    
    private Map<String, Set<String>> generateFNSampleToGenes(Map<String, Set<String>> sampleToGenes,
                                                             double recall,
                                                             RandomDataImpl randomizer) {
        Map<String, Set<String>> rtnSampleToGenes = new HashMap<String, Set<String>>();
        for (String sample : sampleToGenes.keySet()) {
            Set<String> genes = sampleToGenes.get(sample);
            int size = (int)(genes.size() * recall);
            if (size == 0) {
                rtnSampleToGenes.put(sample, new HashSet<String>());
            }
            else {
                Set<String> rtnGenes = MathUtilities.randomSampling(genes, 
                                                                    size,
                                                                    randomizer);
                rtnSampleToGenes.put(sample, rtnGenes);
            }
        }
        return rtnSampleToGenes;
    }
    
    /**
     * Check for false negative caused lower gene number in each sample.
     * Note: this method is copied from another method checkModuleStabilty() for 
     * quick run. This should not be encouraged.
     * @throws Exception
     */
    @Test
    public void checkModuleStabilityForFN() throws Exception {
        // The following statement is used to get the total number of pathways
        AnnotationHelper annotationHelper = new AnnotationHelper();
        Map<String, Set<String>> geneToPathways = annotationHelper.loadProteinNameToTermsMap(AnnotationType.Pathway);
        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
        int totalPathway = pathwayToGenes.size();
        System.out.println("Total pathways: " + totalPathway);
        
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> totalGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> sampleToGenes = loadSampleToMutatedGenes();
        // Used for a random sample set
//        sampleToGenes = generateRandomSamples(sampleToGenes, totalGenes);
        Set<String> selectedGenes = getSelectedGenes(sampleToGenes, totalGenes);
        
        List<Set<String>> clusters = cluster(fis, sampleToGenes);
        Set<String> pathways = annotateClusters(clusters);
        System.out.println("Number of clusters: " + clusters.size());
        System.out.println("Total pathways: " + pathways.size());
        
        // Want to do a sample based random sampling
        RandomDataImpl randomizer = new RandomDataImpl();
        double[] percentages = new double[] {0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20};
        int permutation = 10;
        SummaryStatistics moduleStats = new SummaryStatistics();
        SummaryStatistics sharedModuleStats = new SummaryStatistics();
        SummaryStatistics sharedModulePvalueStats = new SummaryStatistics();
        SummaryStatistics pathwayStats = new SummaryStatistics();
        SummaryStatistics pathwaySharedStats = new SummaryStatistics();
        SummaryStatistics pathwayPValueStats = new SummaryStatistics();
        SummaryStatistics geneStats = new SummaryStatistics();
        StringBuilder builder = new StringBuilder();
        builder.append("Percentage\t" +
                "Selected_Genes\tSD\t" + 
                "Modules\tSD\t" + 
                "Shared_Module\tSD\t" +
                "pValue_Mean\tSD\t" +
                "Annotations\tSD\t" + 
                "Shared_Annotations\tSD\t" +
                "pValue_Mean\tSD\n");
        builder.append("1.00\t");
        builder.append(selectedGenes.size()).append("\t0\t");
        builder.append(clusters.size()).append("\t0\t");
        builder.append("1\t0\t");
        builder.append("0\t0\t");
        builder.append(pathways.size()).append("\t0\t");
        builder.append(pathways.size()).append("\t0\t");
        builder.append("0\t0\n");
        for (double percentage : percentages) {
            System.out.println("\nRecall: " + percentage);
            moduleStats.clear();
            sharedModuleStats.clear();
            sharedModulePvalueStats.clear();
            pathwayStats.clear();
            pathwaySharedStats.clear();
            pathwayPValueStats.clear();
            geneStats.clear();
            for (int i = 0; i < permutation; i++) {
                Map<String, Set<String>> sampleToGenes1 = generateFNSampleToGenes(sampleToGenes,
                                                                                  percentage,
                                                                                  randomizer);
                // Gene sharing without considering clusters
                Set<String> selectedGenes1 = getSelectedGenes(sampleToGenes1, totalGenes);
                geneStats.addValue(selectedGenes1.size());
                // Clusters sharing
                List<Set<String>> clusters1 = cluster(fis, sampleToGenes1);
                moduleStats.addValue(clusters1.size());
                double[] shared = calculateSimilarModules(clusters,
                                                          clusters1,
                                                          totalGenes.size());
                sharedModuleStats.addValue(shared[0]);
                sharedModulePvalueStats.addValue(shared[1]);
                // Annotation sharing
                Set<String> pathways1 = annotateClusters(clusters1);
                double[] pathwayShared = calculateSimilarAnnotations(pathways,
                                                                     pathways1,
                                                                     totalPathway);
                pathwayStats.addValue(pathways1.size());
                pathwaySharedStats.addValue(pathwayShared[0]);
                pathwayPValueStats.addValue(pathwayShared[1]);
            }
            builder.append(percentage + "\t" + 
                    geneStats.getMean() + "\t" + geneStats.getStandardDeviation() + "\t" +
                    moduleStats.getMean() + "\t" + moduleStats.getStandardDeviation() + "\t" +
                    sharedModuleStats.getMean() + "\t" + sharedModuleStats.getStandardDeviation() + "\t" + 
                    sharedModulePvalueStats.getMean() + "\t" + sharedModulePvalueStats.getStandardDeviation() + "\t" + 
                    pathwayStats.getMean() + "\t" + pathwayStats.getStandardDeviation() + "\t" +
                    pathwaySharedStats.getMean() + "\t" + pathwaySharedStats.getStandardDeviation() + "\t" + 
                    pathwayPValueStats.getMean() + "\t" + pathwayPValueStats.getStandardDeviation() + "\n");
        }
        System.out.println(builder.toString());
    }
    
}
