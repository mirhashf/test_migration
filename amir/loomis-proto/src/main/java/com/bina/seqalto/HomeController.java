// Copyright (C) 2012 Bina Technologies Inc.

package com.bina.seqalto;

import java.io.IOException;

import javax.annotation.Resource;
import javax.inject.Inject;
import javax.ws.rs.FormParam;
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;

import org.apache.zookeeper.KeeperException;
import org.apache.zookeeper.KeeperException.Code;
import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.hibernate.SessionFactory;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * Controller bean to respond to JSON API Calls. Built on top of Jersey
 * 
 * @author Amirhossein Kiani (amir@binatechnologies.com)
 */
@Path("/api")
@Component
@Scope("request")
public class HomeController {
  @Resource(name = "sessionFactory")
  @Inject
  public SessionFactory sessionFactory;
  @Inject
  public Leader leader;

  @GET
  @Path("run_alignment")
  @Produces(MediaType.APPLICATION_JSON)
  public Object runAligment() throws Exception {
    JobRequest job = new JobRequest();
    job.setPacketSizeInLines(4);
    job.setProperty(JobParameters.FASTQ_PAIR_1,
        "/Users/amir/bina/analysis/425991/results/samfile.sam");
    job.setProperty(JobParameters.FASTA_FILE, "/Users/amir/bina/genome/hg19/genome_all.fa.fai");
    job.setStreaming(true);
    leader.process(job);
    return job;
  }

  @POST
  @Path("run_alignment_trex")
  @Produces(MediaType.APPLICATION_JSON)
  public Object runAligmentTrex(@FormParam("packet_size") int packetSize,
      @FormParam("aligner_input") String aligner_input, @FormParam("genome") String fastaIndex,
      @FormParam("seqalto_index") String seqaltoIndex) throws Exception {
    if (packetSize == 0) {
      return "Packet size can't be zero!";
    }
    JobRequest job = new JobRequest();
    job.setPacketSizeInLines(packetSize);
    job.setProperty(JobParameters.FASTQ_PAIR_1, aligner_input);
    job.setProperty(JobParameters.FASTA_FILE, SpringPropertiesUtil.getProperty("genome.folder")
        + "/" + fastaIndex);
    job.setProperty(JobParameters.SEQALTO_INDEX, SpringPropertiesUtil.getProperty("genome.folder")
        + "/" + seqaltoIndex);
    job.setStreaming(true);
    leader.process(job);
    return job;
  }

  @GET
  @Path("kill_job/{job_id}")
  @Produces(MediaType.APPLICATION_JSON)
  public Object killJob(@PathParam("job_id") String jobId) throws JsonParseException,
      JsonMappingException, KeeperException, InterruptedException, IOException {
    try {
      leader.killJob(jobId);
    } catch (KeeperException ex) {
      if (ex.code().equals(Code.NONODE)) {
        return "Job " + jobId + " does not exist!";
      }
    }

    return jobId;

  }

  @GET
  @Path("areyouleader")
  @Produces(MediaType.APPLICATION_JSON)
  public boolean areYouLeader() {
    return ZooKeeperService.getIsLeader();
  }


  public void setSessionFactory(SessionFactory sessionFactory) {
    this.sessionFactory = sessionFactory;
  }

  public SessionFactory getSessionFactory() {
    return sessionFactory;
  }
}
