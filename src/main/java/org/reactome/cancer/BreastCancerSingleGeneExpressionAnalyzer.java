/*
 * Created on Jun 30, 2014
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.GeneExpressionDataSet;
import org.reactome.r3.util.InteractionUtilities;

/**
 * This class is used to check single gene expression for breast cancers.
 * @author Gwu
 *
 */
public class BreastCancerSingleGeneExpressionAnalyzer {
//	private final String DIR_NAME = "/Users/Gwu/Desktop/keli/";
	private final String DIR_NAME = "/Users/gwu/Documents/wgm/work/Keli/GSE18229/";
//	private final String DIR_NAME = "/Users/Gwu/Desktop/keli/GSE18229/";
//	private final String DIR_NAME = "/Users/Gwu/Desktop/keli/GSE41748/";
//	private final String DIR_NAME = "/Users/gwu/Documents/wgm/work/Keli/GSE41748/";
	private GSESoftTextDataHandler dataHandler;
	private FileUtility fu;

	public BreastCancerSingleGeneExpressionAnalyzer() {
		dataHandler = new GSESoftTextDataHandler();
		fu = new FileUtility();
	}
	
	/**
	 * The output from this file is saved as ProbeIdAnnotations.txt after some cleaning up.
	 * @throws IOException
	 */
	@Test
	public void checkGSE18229Platforms() throws IOException {
		// Use a pre-processed platform files copied from the original family file.
		String fileName = DIR_NAME + "GSE18229_family_platforms.soft";
		boolean isMouse = false;
//		String fileName = DIR_NAME + "GSE41748_family.soft";
//		boolean isMouse = true;
		List<String> platforms = dataHandler.loadPlatforms(fileName);
		System.out.println("Total platforms: " + platforms.size());
		for (String platform : platforms)
			System.out.println(platform);
		fu.setInput(fileName);
		String line = null;
		List<String> genes = getTargetGenes();
		while ((line = fu.readLine()) != null) {
			if (line.startsWith("^PLATFORM ="))
				System.out.println(line);
			else if (line.startsWith("!platform_table_begin")) {
				line = fu.readLine();
				System.out.println(line);
//				line = fu.readLine();
//				System.out.println(line);
			}
			else {
				if (isMouse)
					line = line.toUpperCase();
				// This should be more reliable by checking individual tokens
                String[] tokens = line.split("\t");
				for (String gene : genes) {
				    for (String token : tokens) {
				        if (token.equals(gene)) {
				            System.out.println(line);
				            break;
				        }
				    }
//					if (line.contains(gene)) {
//						System.out.println(line);
//						break;
//					}
				}
			}
		}
		fu.close();
//		Map<String, String> probeIdToGene = dataHandler.loadProbeIdToGene(fileName);
//        System.out.println("Total probe ids: " + probeIdToGene.size());
	}

    private List<String> getTargetGenes() {
        // List of genes to be pulled out
		List<String> genes = new ArrayList<String>();
		genes.add("MFNG");
		genes.add("NOTCH4");
		genes.add("PIK3CG");
		genes.add("MYC");
        return genes;
    }
    
    private List<String> getMouseTargetGenes() {
        List<String> genes = new ArrayList<String>();
        genes.add("Mfng");
        genes.add("Notch4");
        genes.add("Pik3cg");
        return genes;
    }
    
	@Test
	public void processCancerTypesFromGSE18229() throws IOException {
		String fileName = DIR_NAME + "GSE18229_family.soft";
		String featureName = "!Sample_characteristics_ch2 = pam50 predictions plus claudin-low classification (cell line predictor)";
		Map<String, String> sampleToFeature = dataHandler.loadSampleFeature(fileName, 
																	        featureName, 
																	        ":");
		System.out.println("sampleToFeature: " + sampleToFeature.size());
		for (String sample : sampleToFeature.keySet()) {
			System.out.println(sample + "\t" + sampleToFeature.get(sample));
		}
	}
	
	@Test
	public void processCancerTypesFromGSE41748() throws IOException {
		String fileName = DIR_NAME + "GSE41748_family.soft";
		String featureName = "!Sample_characteristics_ch1 = pathology";
		Map<String, String> sampleToFeature = dataHandler.loadSampleFeature(fileName, featureName, ":");
		System.out.println("sampleToFeature: " + sampleToFeature.size());
		for (String sample : sampleToFeature.keySet())
			System.out.println(sample + "\t" + sampleToFeature.get(sample));
	}
	
	   @Test
	    public void processGSE41748InProbes() throws IOException {
	        List<String> targetGenes = getMouseTargetGenes();
	        String fileName = DIR_NAME + "GSE41748_family.soft";
	        List<String> escaped = new ArrayList<String>();
	        escaped.add("GPL16199");
	        Map<String, String> platProbGenes = dataHandler.loadProbeIdToGene(fileName, escaped);
	        List<String> selectedProbes = new ArrayList<String>();
	        for (String key : platProbGenes.keySet()) {
	            String gene = platProbGenes.get(key);
	            if (targetGenes.contains(gene)) {
	                System.out.println(key + "\t" + gene);
	                selectedProbes.add(key);
	            }
	        }
	        Map<String, List<String>> platformToProbeIds = new HashMap<String, List<String>>();
	        for (String selectedProbe : selectedProbes) {
	            int index = selectedProbe.indexOf(".");
	            String platform = selectedProbe.substring(0, index);
	            String id = selectedProbe.substring(index + 1);
	            List<String> ids = platformToProbeIds.get(platform);
	            if (ids == null) {
	                ids = new ArrayList<String>();
	                platformToProbeIds.put(platform, ids);
	            }
	            ids.add(id);
	        }
	        File[] files = new File(DIR_NAME).listFiles();
	        List<String> samples = new ArrayList<String>();
	        GSEMatrixDataHandler dataHandler = new GSEMatrixDataHandler();
	        for (File file : files) {
	            String name = file.getName();
	            if (!name.endsWith("_matrix.txt") || name.contains("GPL16199")) 
	                continue;
	            //              if (!name.equals("GSE18229-GPL5325_series_matrix.txt"))
	            //                  continue;
	            System.out.println("File: " + file.getName());
	            String platform = parsePlatform(name);
	            System.out.println("Platform: " + platform);
	            GeneExpressionDataSet dataSet = dataHandler.loadGeneExpressionDataSet(file.getAbsolutePath());
	            dataSet = dataSet.selectFeatures(platformToProbeIds.get(platform));
	            dataSet.export(DIR_NAME + platform + "_selected_probes.txt");
	        }
	    }
	
	@Test
	public void processGSE41748() throws IOException {
	    List<String> targetGenes = getMouseTargetGenes();
	    String fileName = DIR_NAME + "GSE41748_family.soft";
	    List<String> escaped = new ArrayList<String>();
	    escaped.add("GPL16199");
	    Map<String, String> platProbGenes = dataHandler.loadProbeIdToGene(fileName, escaped);
	    for (String key : platProbGenes.keySet()) {
	        String gene = platProbGenes.get(key);
	        if (targetGenes.contains(gene)) {
	            System.out.println(key + "\t" + gene);
	        }
	    }
	    Map<String, Map<String, String>> platformToProbeIdToGene = splitPlatform(platProbGenes);
	    File[] files = new File(DIR_NAME).listFiles();
	    List<String> samples = new ArrayList<String>();
	    GSEMatrixDataHandler dataHandler = new GSEMatrixDataHandler();
	    List<GeneExpressionDataSet> datasets = new ArrayList<GeneExpressionDataSet>();
	    for (File file : files) {
	        String name = file.getName();
	        if (!name.endsWith("_matrix.txt") || name.contains("GPL16199")) 
	            continue;
	        //	            if (!name.equals("GSE18229-GPL5325_series_matrix.txt"))
	        //	                continue;
	        System.out.println("File: " + file.getName());
	        String platform = parsePlatform(name);
	        System.out.println("Platform: " + platform);
	        GeneExpressionDataSet dataSet = dataHandler.loadGeneExpressionDataSet(file.getAbsolutePath());
	        dataHandler.mapProbesetToGenes(dataSet, platformToProbeIdToGene.get(platform));
	        dataSet = dataHandler.averageValuesForSameGenes(dataSet);
	        dataSet.sampleWiseRankTransformation();
	        dataSet = dataSet.selectFeatures(targetGenes);
	        datasets.add(dataSet);
	    }
	    GeneExpressionDataSet merged = this.dataHandler.mergeDatasetsWithDiffSamples(datasets);
	    merged.export(DIR_NAME + "Merged_Expression_Values_Rank.txt");
	}
	
	private Map<String, Map<String, String>> splitPlatform(Map<String, String> platProbGenes) {
	    Map<String, Map<String, String>> rtn = new HashMap<String, Map<String,String>>();
	    for (String platProb : platProbGenes.keySet()) {
	        String gene = platProbGenes.get(platProb);
	        int index = platProb.indexOf(".");
	        String platform = platProb.substring(0, index);
	        String probeId = platProb.substring(index + 1);
	        Map<String, String> probeIdToGene = rtn.get(platform);
	        if (probeIdToGene == null) {
	            probeIdToGene = new HashMap<String, String>();
	            rtn.put(platform, probeIdToGene);
	        }
	        probeIdToGene.put(probeId, gene);
	    }
	    return rtn;
	}
	
	@Test
	public void grepExpressionValuesForSingleGenesInGSE18229() throws IOException {
	    Map<String, Map<String, String>> platProbGenes = loadPlatformAnnotFile();
	    for (String plat : platProbGenes.keySet()) {
	        Map<String, String> probToGene = platProbGenes.get(plat);
	        System.out.println(plat + ": " + probToGene);
	    }
	    File[] files = new File(DIR_NAME).listFiles();
        List<String> samples = new ArrayList<String>();
        GSEMatrixDataHandler dataHandler = new GSEMatrixDataHandler();
        List<GeneExpressionDataSet> datasets = new ArrayList<GeneExpressionDataSet>();
        for (File file : files) {
            String name = file.getName();
            if (!name.endsWith("_matrix.txt") || name.contains("GPL16199")) 
                continue;
//            if (!name.equals("GSE18229-GPL5325_series_matrix.txt"))
//                continue;
            System.out.println("File: " + file.getName());
            String platform = parsePlatform(name);
            System.out.println("Platform: " + platform);
            GeneExpressionDataSet dataSet = dataHandler.loadGeneExpressionDataSet(file.getAbsolutePath());
            dataHandler.mapProbesetToGenes(dataSet, platProbGenes.get(platform));
            dataSet = dataHandler.averageValuesForSameGenes(dataSet);
            datasets.add(dataSet);
        }
        GeneExpressionDataSet merged = this.dataHandler.mergeDatasetsWithDiffSamples(datasets);
        merged.export(DIR_NAME + "Merged_Expression_Values_030416.txt");
	}
	
	@Test
	public void generateMergedExpressionFile() throws IOException {
	    Map<String, String> sampleToType = loadSampleToSubtype();
//	    String fileName = DIR_NAME + "Merged_Expression_Values.txt";
//	    String fileName = DIR_NAME + "Merged_Expression_Values_Rank.txt";
	    String fileName = DIR_NAME + "GPL2884_selected_probes.txt";
	    List<String[]> list = new ArrayList<String[]>();
	    fu.setInput(fileName);
	    String line = null;
	    int tokenSize = 0;
	    while ((line = fu.readLine()) != null) {
	        String[] tokens = line.split("\t");
	        list.add(tokens);
	        if (tokenSize == 0)
	            tokenSize = tokens.length;
	    }
	    fu.close();
//	    fu.setOutput(DIR_NAME + "Merged_Expression_Values_Types.txt");
//	    fu.setOutput(DIR_NAME + "Merged_Expression_Values_Rank_Types.txt");
	    fu.setOutput(DIR_NAME + "GPL2884_selected_probes_Types.txt");
	    StringBuilder builder = new StringBuilder();
	    for (int i = 0; i < tokenSize; i++) {
	        String sample = null;
	        for (int j = 0; j < list.size(); j++) {
	            String[] tokens = list.get(j);
	            builder.append(tokens[i]).append("\t");
	            if (j == 0)
	                sample = tokens[i];
	        }
	        boolean escape = false;
	        if (i == 0)
	            builder.append("Type");
	        else {
	            String type = sampleToType.get(sample);
	            if (type != null && !type.equals("NA"))
	                builder.append(type);
	            else {
	                //builder.append("");
	                escape = true;
	            }
	        }
	        if (!escape)
	            fu.printLine(builder.toString());;
	        builder.setLength(0);
	    }
	    fu.close();
	}
	
	private Map<String, String> loadSampleToSubtype() throws IOException {
	    String fileName = DIR_NAME + "Sample_List_WithTypes.txt";
	    return fu.loadInteractionPairs(fileName, "\t");
	}
	
	private String parsePlatform(String fileName) {
	    int index = fileName.indexOf("-");
	    int index1 = fileName.indexOf("_");
	    return fileName.substring(index + 1, index1);
	}
	
	private Map<String, Map<String, String>> loadPlatformAnnotFile() throws IOException {
	    String fileName = DIR_NAME + "ProbeIdAnnotations_030416.txt";
	    fu.setInput(fileName);
	    String line = null;
	    String platform = null;
	    int index = 0;
	    List<String> genes = getTargetGenes();
	    Map<String, Map<String, String>> rtn = new HashMap<String, Map<String,String>>();
	    Map<String, String> probeIdToGene = null;
	    while ((line = fu.readLine()) != null) {
	        if (line.startsWith("^PLATFORM")) {
	            index = line.indexOf("=");
	            platform = line.substring(index + 1).trim();
	            probeIdToGene = rtn.get(platform);
	            if (probeIdToGene == null) {
	                probeIdToGene = new HashMap<String, String>();
	                rtn.put(platform, probeIdToGene);
	            }
	        }
	        else if (line.startsWith("ID"))
	            continue;
	        else {
	            String[] tokens = line.split("\t");
	            // Find gene
	            for (String token : tokens) {
	                if (genes.contains(token)) {
	                    probeIdToGene.put(tokens[0], token);
	                    break;
	                }
	            }
	        }
	    }
	    fu.close();
	    return rtn;
	}
	
	@Test
	public void checkSamplesInFilesForGSE18229() throws IOException {
	    File[] files = new File(DIR_NAME).listFiles();
	    List<String> samples = new ArrayList<String>();
	    for (File file : files) {
	        String name = file.getName();
	        if (!name.endsWith("_matrix.txt")) 
	            continue;
	        System.out.println("File: " + file.getName());
	        fu.setInput(file.getAbsolutePath());
	        String line = null;
	        while ((line = fu.readLine()) != null) {
	            if (line.startsWith("\"ID_REF\"")) {
	                System.out.println(line);
	                String[] tokens = line.split("\t");
	                for (int i = 1; i < tokens.length; i++)
	                    samples.add(tokens[i]);
	                break;
	            }
	        }
	        System.out.println("Size: " + samples.size());
	        fu.close();
	    }
	    System.out.println("Total samples: " + samples.size());
	    System.out.println("    In set: " + new HashSet<String>(samples).size());
	}
	
	@Test
	public void checkSamplesInFilesForGSE41748() throws IOException {
	    File[] files = new File(DIR_NAME).listFiles();
	    List<String> samples = new ArrayList<String>();
	    for (File file : files) {
	        String name = file.getName();
	        if (!name.endsWith("_matrix.txt") || name.contains("GPL16199")) // 16199 is a miRNA array
	            continue;
	        System.out.println("File: " + file.getName());
	        fu.setInput(file.getAbsolutePath());
	        String line = null;
	        while ((line = fu.readLine()) != null) {
	            if (line.startsWith("\"ID_REF\"")) {
	                System.out.println(line);
	                String[] tokens = line.split("\t");
	                for (int i = 1; i < tokens.length; i++)
	                    samples.add(tokens[i]);
	                break;
	            }
	        }
	        System.out.println("Size: " + samples.size());
	        fu.close();
	    }
	    System.out.println("Total samples: " + samples.size());
	    String fileName = DIR_NAME + "GSE41748_family.soft";
		String featureName = "!Sample_title";
		Map<String, String> sampleToTitle = dataHandler.loadSampleFeature(fileName, 
				featureName, 
				"=");
		Map<String, String> sampleToType = fu.loadInteractionPairs(DIR_NAME + "SampleToType.txt", "\t");
		Set<String> titles = new HashSet<String>();
		for (String sample : samples) {
			sample = InteractionUtilities.removeQuotationMarks(sample);
			String title = sampleToTitle.get(sample);
			String type = sampleToType.get(sample);
			System.out.println(sample + "\t" + title + "\t" + type);
			titles.add(title);
		}
		System.out.println("Total titles: " + titles.size());
	}
	
	@Test
	public void checkGSE18229Samples() throws IOException {
		String fileName = DIR_NAME + "Sample_List.txt";
		List<String> samples1 = new ArrayList<String>();
		List<String> samples2 = new ArrayList<String>();
		fu.setInput(fileName);
		String line = null;
		while ((line = fu.readLine()) != null) {
			String[] tokens = line.split("\t");
			samples1.add(tokens[0]);
			if (tokens.length == 4)
				samples2.add(tokens[2]);
		}
		fu.close();
		System.out.println("List 1: " + samples1.size());
		System.out.println("List 2: " + samples2.size());
		Set<String> notShared = new HashSet<String>(samples1);
		notShared.removeAll(samples2);
		System.out.println("Not shared: " + notShared.size());
		List<String> list = new ArrayList<String>(notShared);
		Collections.sort(list);
		for (String sample : list)
			System.out.println(sample);
	}
	
	@Test
	public void checkSampleTitles() throws IOException {
//		String fileName = DIR_NAME + "GSE18229_family.soft";
		String fileName = DIR_NAME + "GSE41748_family.soft";
		String featureName = "!Sample_title";
		Map<String, String> sampleToFeature = dataHandler.loadSampleFeature(fileName, 
																	        featureName, 
																	        "=");
		System.out.println("sampleToFeature: " + sampleToFeature.size());
		System.out.println("Total samples: " + sampleToFeature.size());
		System.out.println("Sample titles: " + new HashSet<String>(sampleToFeature.values()).size());
		for (String sample : sampleToFeature.keySet()) {
			System.out.println(sample + "\t" + sampleToFeature.get(sample));
		}
	}
	
	@Test
	public void processGSE42460SuppTable() throws IOException {
	    String dir = "/Users/gwu/Documents/wgm/work/Keli/GSE42460/";
	    File[] files = new File(dir).listFiles();
	    int totalType = 0;
	    List<String> mouseTargetGenes = getMouseTargetGenes();
	    StringBuilder builder = new StringBuilder();
	    builder.append("Type\tGeneId\tGeneName\tSAMFoldChange\tSAM_q_value(%)");
	    System.out.println(builder.toString());
	    for (File file : files) {
	        String fileName = file.getName();
	        if (!fileName.endsWith(".txt")) {
	            continue;
	        }
	        int index = fileName.indexOf(".");
	        String type = fileName.substring(0, index);
//	        System.out.println(type);
	        totalType ++;
	        fu.setInput(file.getAbsolutePath());
	        String line = fu.readLine();
	        while ((line = fu.readLine()) != null) {
	            String[] tokens = line.split("\t");
	            if (mouseTargetGenes.contains(tokens[1])) {
	                builder.setLength(0);
	                builder.append(type);
	                for (int i = 0; i < 4; i++)
	                    builder.append("\t").append(tokens[i]);
	                System.out.println(builder.toString());
	            }
	        }
	        fu.close();
	    }
	    System.out.println("Total types: " + totalType);
	}
}
