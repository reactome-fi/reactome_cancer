/*
 * Created on Jun 3, 2009
 *
 */
package org.reactome.cancer;

/**
 * This class is used to model CNV data set.
 * @author wgm
 *
 */
public class CopyNumberVariation {
    
    private String gene;
    private int value; // 0, 1, 2, -1, -2
    
    public CopyNumberVariation() {
    }
    
    public String getGene() {
        return gene;
    }

    public void setGene(String gene) {
        this.gene = gene;
    }

    public int getValue() {
        return value;
    }

    public void setValue(int value) {
        this.value = value;
    }
}
