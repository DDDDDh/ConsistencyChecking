package cn.edu.nju.moon.consistency.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import cn.edu.nju.moon.consistency.model.operation.GenericOperation;

/**
 * @author hengxin
 * @date 2012-12-6
 * 
 * @description construct a raw observation (RawObservation) randomly
 */
public class RandomRawObservationConstructor
{
	// number of processes (RawProcess)
	private int processNum = 5;
	// number of variables involved (from 'a' to 'a' + (variableNum - 1))
	private int variableNum = 5;
	// domain of all variables ([0, valueRange))
	private int valueRange = 10;
	// number of operations in all (Operation)
	private int opNum = 30;


	/**
	 * Default value:
	 * variableNum = 5;
	 * valueRange = [0,10);
	 * opNum = 30; (number of processes)
	 * processNum = 5;
	 */
	public RandomRawObservationConstructor(int processNum, int variableNum, int valueRange, int opNum)
	{
		this.processNum = processNum;
		this.variableNum = variableNum;
		this.valueRange = valueRange;
		this.opNum = opNum;
	}
	
	public void construct()
	{
		
	}
	
	/**
	 * @return List of Operation.
	 *   Constraint: Write distinct values for a variable.
	 */
	private List<GenericOperation> constructOperationList()
	{
		List<GenericOperation> opList = new ArrayList<GenericOperation>();
		Random random = new Random();
		GenericOperation op = null;
		int loop = 0;

		for(int i=0;i<this.opNum;)
		{
			loop++;

			op = this.consturctOperation();

			// if it is a Read operation, generate a corresponding Write operation
			if(op.isReadOp())
			{
				opList.add(op);
				i++;

				GenericOperation wop = new GenericOperation(GlobalData.WRITE, op.getVariable(), op.getValue());

				if(i == this.opNum)
				{
					if(! opList.contains(wop))
					{
						opList.remove(op);
						i--;
					}
				}
				else
				{
					if(! opList.contains(wop))
					{
						opList.add(wop);
						i++;
					}
				}
			}
			else // Write Operation
			{
				if(! opList.contains(op) && random.nextInt(3) == 0)
				{
					opList.add(op);
					i++;
				}
			}
		}

//		System.out.println("#" + loop + " loop for " + this.opNum + " Operations");

		Collections.shuffle(opList);

		for(GenericOperation gop : opList)
			if(! gop.isReadOp())
				GlobalData.WRITEPOOL.put(gop.toString(), gop);

		return opList;

	}

	/**
	 * construct generic operation randomly
	 * 
	 * @return {@link GenericOperation} object
	 */
	private GenericOperation consturctOperation()
	{
		Random typeRandom = new Random();
		Random variableRandom = new Random();
		Random valueRandom = new Random();

		// generate a random operation
		int type = 0;
		if(typeRandom.nextBoolean())
			type = GlobalData.READ;
		else
			type = GlobalData.WRITE;

		String var = String.valueOf((char) ('a' + variableRandom.nextInt(this.variableNum)));
		int val = valueRandom.nextInt(this.valueRange);

		return new GenericOperation(type, var, val);
	}
}