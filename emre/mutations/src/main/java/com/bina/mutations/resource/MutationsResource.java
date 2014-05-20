package com.bina.mutations.resource;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.ws.rs.DefaultValue;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.QueryParam;
import javax.ws.rs.core.MediaType;

import com.bina.mutations.model.Annotation;
import com.bina.mutations.model.Annotation.Gene;
import com.bina.mutations.model.Annotation.Impact;
import com.bina.mutations.model.Mutation;
import com.bina.mutations.model.QueryResult;

@Path("/mutations")
@Produces(MediaType.APPLICATION_JSON)
public class MutationsResource {

  private static final Random random = new Random();
  private static final Map<Integer, List<Annotation>> annMap = new HashMap<>();
  private static final Mutation[] snps = generateRandomSnps();

  private static Mutation[] generateRandomSnps() {
    int numMutations = 100;
    int pos = 10;
    Mutation[] result = new Mutation[numMutations];
    for (int i = 1; i <= numMutations; i++) {
      Mutation mut = new Mutation();
      mut.setId(i);
      pos += 10;
      mut.setPosition(pos);
      mut.setChromosome(Mutation.Chromosome.values()[random.nextInt(Mutation.Chromosome.values().length)]);
      int refOrd = random.nextInt(Mutation.Base.values().length);
      mut.setReference(Mutation.Base.values()[refOrd].name());
      int altOrd = random.nextInt(Mutation.Base.values().length);
      if (altOrd == refOrd) {
        if (refOrd == Mutation.Base.values().length - 1) {
          altOrd = 0;
        } else {
          altOrd = refOrd + 1;
        }
      }
      mut.setAlternate(Mutation.Base.values()[altOrd].name());
      mut.setReadDepth(10 + random.nextInt(60));

      int numAnnotations = 1 + random.nextInt(5);
      List<Annotation> annotations = generateRandomAnnotations(numAnnotations);
      annMap.put(mut.getId(), annotations);

      result[i - 1] = mut;
    }
    return result;
  }

  private static List<Annotation> generateRandomAnnotations(int numAnnotations) {
    List<Annotation> annotations = new ArrayList<>(numAnnotations);
    for (int x = 0; x < numAnnotations; x++) {
      Annotation ann = new Annotation();
      ann.setImpact(Impact.values()[random.nextInt(Impact.values().length)]);
      ann.setGene(Gene.values()[random.nextInt(Gene.values().length)]);
      ann.setTranscriptId("ENST" + String.valueOf(1000 + random.nextInt(1000)));
      annotations.add(ann);
    }
    return annotations;
  }

  @GET
  public QueryResult getMutations(@DefaultValue("1") @QueryParam("page") int page,
                                  @DefaultValue("10") @QueryParam("pagesize") int pageSize) {
    return new QueryResult(Arrays.asList(Arrays.copyOfRange(snps, (page - 1) * pageSize, page * pageSize)),
                           page,
                           pageSize,
                           snps.length);
  }

  @GET
  @Path("{id}/annotations")
  public List<Annotation> getAnnotation(@PathParam("id") Integer id) {
    return annMap.get(id);
  }
}
