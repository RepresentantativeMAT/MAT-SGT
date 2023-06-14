/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package br.ufsc.model;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.Objects;

/**
 *
 * @author vanes
 */
public class TemporalAspect {
    
    private Date startTime;
    private Date endTime;
    SimpleDateFormat formatDate = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
    SimpleDateFormat formatTime = new SimpleDateFormat("HH:mm:ss");
    private boolean dailyInfo;

    public TemporalAspect(){}
     
    public TemporalAspect(Date startTime) {
        this.startTime = startTime;
    }
    

    public TemporalAspect(Date startTime, Date endTime) {
        this.startTime = startTime;
        this.endTime = endTime;
    }
    
    public TemporalAspect(int startTimeMinutes) {
        this.startTime = Util.convertMinutesToDate(startTimeMinutes);
        setDailyInfo(true);
    }

    public Date getStartTime() {
        return startTime;
    }

    public void setStartTime(Date startTime) {
        this.startTime = startTime;
    }
    public void setStartTime(int startTimeMinutes) {
        this.startTime = Util.convertMinutesToDate(startTimeMinutes);
    }

    public Date getEndTime() {
        return endTime;
    }

    public void setEndTime(Date endTime) {
        this.endTime = endTime;
    }
    
    public void setEndTime(int endTimeMinutes) {
        this.endTime = Util.convertMinutesToDate(endTimeMinutes);
    }

    @Override
    public String toString() {
        return formatDate.format(startTime) +" - "+ (endTime != null? formatDate.format(endTime):"");
    }

    /**
     * Analyse if the time is contained in this interval
     * @param time
     * @return 
     */
    public boolean isInInterval(Date time){
        Calendar start = Calendar.getInstance();
             start.set(Calendar.HOUR_OF_DAY, startTime.getHours());
             start.set(Calendar.MINUTE, startTime.getMinutes());
             
             Calendar compared = Calendar.getInstance();
             compared.set(Calendar.HOUR_OF_DAY, time.getHours());
             compared.set(Calendar.MINUTE, time.getMinutes());
             
             
        try{
            getEndTime();
            
            // To analyse if the current point is in the interval -- print current point x intervals             
             Calendar end = Calendar.getInstance();
             end.set(Calendar.HOUR_OF_DAY, endTime.getHours());
             end.set(Calendar.MINUTE, endTime.getMinutes());
             
             
            return (time.compareTo(startTime) >0 && time.compareTo(endTime)<0 
//                    || time.compareTo(startTime) == 0 
//                    || time.compareTo(endTime) ==0 );
                    || start.compareTo(compared) == 0
                    || end.compareTo(compared) == 0);
            
        }catch (NullPointerException e){
            return start.compareTo(compared) == 0;
        }
    }
    
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Objects.hashCode(this.startTime);
        hash = 29 * hash + Objects.hashCode(this.endTime);
        hash = 29 * hash + Objects.hashCode(this.formatDate);
        hash = 29 * hash + Objects.hashCode(this.formatTime);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final TemporalAspect other = (TemporalAspect) obj;
        if (!Objects.equals(this.startTime, other.startTime)) {
            return false;
        }
        if (!Objects.equals(this.endTime, other.endTime)) {
            return false;
        }
        if (!Objects.equals(this.formatDate, other.formatDate)) {
            return false;
        }
        if (!Objects.equals(this.formatTime, other.formatTime)) {
            return false;
        }
        return true;
    }

    public boolean isDailyInfo() {
        return dailyInfo;
    }

    public void setDailyInfo(boolean dailyInfo) {
        if(dailyInfo == true)
            formatDate = new SimpleDateFormat("HH:mm");
        this.dailyInfo = dailyInfo;
    }
    
    
    
}
