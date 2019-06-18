/*
 * Created on May 11, 2010
 *
 */
package org.reactome.cancer;

/**
 * This simple class is used to model clinical information from TCGA.
 * @author wgm
 *
 */
public class TCGAClinicalInfo {
    private String sample;
    private Boolean osEvent;
    private Double osDuration;
    
    public TCGAClinicalInfo() {
        
    }
    
    public void setSample(String sample) {
        this.sample = sample;
    }
    
    public String getSample() {
        return this.sample;
    }

    public Boolean getOsEvent() {
        return osEvent;
    }

    public void setOsEvent(Boolean osEvent) {
        this.osEvent = osEvent;
    }

    public Double getOsDuration() {
        return osDuration;
    }

    public void setOsDuration(Double osDuration) {
        this.osDuration = osDuration;
    }
    
    
    
}
