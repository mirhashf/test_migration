/**
 * Roman Number Adder Test
 * @author Kuang Sheng
 * @summary this program tests all the defined methods in romanNumberAdder class
 */

import static org.junit.Assert.*;

import org.junit.Test;

public class romanNumberAdderTest {
	
	@Test
	public void test_add(){
		romanNumberAdder rna = new romanNumberAdder();
		assertEquals("IX",rna.add("","IX"));
		assertEquals("XII",rna.add("III","IX"));
		assertEquals("CCX",rna.add("CLI","LIX"));
		assertEquals("XX",rna.add("XI","VIIII"));
		assertEquals("MCCCXIII",rna.add("MCLXXIV","CXXXIX"));
	}
	
	@Test
	public void test_convertSubtractiveToAdditive(){
		romanNumberAdder rna = new romanNumberAdder();
		assertEquals("",rna.convertSubtractiveToAdditive(""));
		assertEquals("DCCCCCCCCLXXXXXXXXVIIIIIIII",rna.convertSubtractiveToAdditive("CMCDXCXLIXIV"));
	}
	
	@Test
	public void test_combineAdditive(){
		romanNumberAdder rna = new romanNumberAdder();
		assertEquals("",rna.combineAdditive(""));
		assertEquals("MCCCXIII",rna.combineAdditive("MCCLXXXXXVIIIIIIII"));
		assertEquals("MCX",rna.combineAdditive("CCCCCCCCCCXXXXXXXXXXIIIIIIIIII"));
		assertEquals("MDCLXXV",rna.combineAdditive("DDDLLLVVVVV"));
	}
	
	@Test
	public void test_characterArrayToString(){
		romanNumberAdder rna = new romanNumberAdder();
		assertEquals("",rna.characterArrayToString(new Character[]{}));
		assertEquals("MCCCXIII",rna.characterArrayToString(
				    new Character[]{'M','C','C','C','X','I','I','I'}));
	}
	
	
	

}
