/*
 * Created on Apr 13, 2017
 *
 */
package org.reactome.r3;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;

//import jonelo.jacksum.JacksumAPI;
//import jonelo.jacksum.algorithm.AbstractChecksum;

/**
 * Generate a map between different databases based on different criteria.
 * @author gwu
 *
 */
public class ChemicalAnalyzer {
    private FileUtility fu = new FileUtility();
    
    /**
     * Default constructor.
     */
    public ChemicalAnalyzer() {
    }
    
    private Set<String> getReactomeChEBIs() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_59_plus_i",
                                            "root",
                                            "macmysql01");
        Set<String> ids = new HashSet<>();
        Collection<GKInstance> c = dba.fetchInstancesByClass(ReactomeJavaConstants.ReferenceMolecule);
        Collection<GKInstance> all = new HashSet<>(c);
        c = dba.fetchInstancesByClass(ReactomeJavaConstants.ReferenceGroup);
        all.addAll(c);
        for (GKInstance inst : all) {
            GKInstance dbName = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            if (dbName == null)
                continue;
            if (!dbName.getDisplayName().equals("ChEBI"))
                continue;
            String identifier = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
            ids.add(identifier);
        }
        return ids;
    }
    
    /**
     * This map is based on InChI information.
     * @throws Exception
     */
    @Test
    public void mapChEBIAndPDB() throws Exception {
        String dirName = "datasets/chemicals/";
        
        Map<String, String> inChiToChebi = new HashMap<>();
        String chebiFileName = dirName + "chebi/chebiId_inchi.tsv";
        fu.setInput(chebiFileName);
        String line = fu.readLine(); // There is a header
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            inChiToChebi.put(tokens[1], tokens[0]);
        }
        fu.close();
        System.out.println("Total InChi in ChEBI: " + inChiToChebi.size());
        
        String pdbFileName = dirName + "pdb/Components-inchi.ich.txt";
        fu.setInput(pdbFileName);
        // There is no header in this file
        Map<String, String> inChiToPDB = new HashMap<>();
        while ((line = fu.readLine()) != null) {
            // Some format error in the downloaded file
//            System.out.println(line);
            if (line.length() == 0)
                continue;
            String[] tokens = line.split("\t");
            if (tokens.length < 2)
                continue;
            inChiToPDB.put(tokens[0], tokens[1]);
        }
        System.out.println("Total InChi in PDB: " + inChiToPDB.size());
        
        // Map PDB chemicals to ChEBI ids
        int count = 0;
        String outFileName = dirName + "PDBChemicalsToChEBIIds_041317.txt";
        fu.setOutput(outFileName);
        fu.printLine("PDB_Chemicals\tChEBI_ID");
        Set<String> mappedChEBI = new HashSet<>();
        for (String pdbInChi : inChiToPDB.keySet()) {
            String ebiId = inChiToChebi.get(pdbInChi);
            if (ebiId == null) {
                count ++;
                continue;
            }
            fu.printLine(inChiToPDB.get(pdbInChi) + "\t" + ebiId);
            mappedChEBI.add(ebiId);
        }
        fu.close();
//        System.out.println("Total InChi in PDB: " + inChiToPDB.size());
        System.out.println("Total not mapped: " + count);
        Set<String> chebiInReactome = getReactomeChEBIs();
        System.out.println("Total ChEBI in reactome: " + chebiInReactome.size());
        chebiInReactome.retainAll(mappedChEBI);
        System.out.println("\tMapped to PDB chemicals: " + chebiInReactome.size());
    }
    
    
    
    /**
     * The following cannot work. Need to figure out a better way.
     * @throws Exception
     */
    @Test
    public void testInChiKey() throws Exception {
//        AbstractChecksum checksum = JacksumAPI.getChecksumInstance("sha-256");
//        checksum.setEncoding(AbstractChecksum.DEC);
//        String InChI = "1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-5,13,15-20H,6H2/t13-,15-/m1/s1";
//        System.out.println(InChI);
//        checksum.reset();
//        checksum.update(InChI.getBytes());
//        System.out.println("InChIKey: " + checksum.getFormattedValue());
    }
    
}
