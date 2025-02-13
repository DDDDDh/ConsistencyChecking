package src.cn.edu.nju.moon.consistency.schedule;

import java.util.Arrays;

import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;

/**
 * @description 
 * In <bf>weak</bf> consistency models, different Processes may hold different {@link View}s
 * satisfying consistency condition.
 *
 * A {@link WeakSchedule} is a group of {@link View}s seen by different Processes.
 * 
 * @author hengxin
 * @date 2013-1-9
 */
public class WeakSchedule implements ISchedule
{

	private int procNum = -1;	/** number of processes */
	private View[] views = null;	/** one {@link View} for each Process */
	private Boolean[] binaryViews = null;
	/**
	 * Constructor: allocate one view for each process 
	 * @param procNum number of processes
	 */
	public WeakSchedule(int procNum)
	{
		this.procNum = procNum;
		this.views = new View[procNum];
		this.binaryViews = new Boolean[procNum];
	}
	
	/**
	 * check whether the {@link WeakSchedule} (i.e., {@link #views}) satisfies PRAM Consistency.
	 * {@link #views} satisfies PRAM Consistency if and only if each {@link View} satisfies 
	 * PRAM Consistency
	 * 
	 * @return true, if {@link #views} satisfies PRAM Consistency; false, otherwise.
	 */
	@Override
	public boolean valid()
	{
		for (View v : this.views)
		{	
			if (v == null)
				return false;
			if (! v.self_check())
				return false;
		}
		return true;
	}

	/**
	 * construct a possible {@link View} from @param bob
	 * @param bob {@link BasicObservation} from which a possible {@link View} is constructed
	 */
	public void constructView(BasicObservation bob)
	{
		this.views[bob.getMasterPid()] = new View(bob);
	}

	/**
	 * get the binary schedule in the sense that the concrete {@link View}
	 * is not significant, but its existence is.
	 * 
	 * @return binary schedule: 1 for consistent; 0 for inconsistent.
	 */
	public Boolean[] constructBinarySchedule()
	{
		for (int index = 0; index < this.procNum; index++)
			if (this.views[index] == null)
				this.binaryViews[index] = false;
			else
				this.binaryViews[index] = true;
		
		return this.binaryViews;
	}
	
	/**
	 * compare two {@link WeakSchedule}s to see whether they are compatible to each other
	 * @param ws another {@link WeakSchedule} to compare with
	 * @return true, if the two {@link WeakSchedule} are compatible to each other; false, otherwise.
	 */
	@Override
	public boolean compare(ISchedule ws)
	{
		if ( ! (ws instanceof WeakSchedule))
			return false;
		
		return Arrays.equals(this.constructBinarySchedule(), ((WeakSchedule) ws).constructBinarySchedule());
	}
	
	/**
	 * @return String form of {@link #views}: one line for each view 
	 */
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		for (int index = 0; index < this.procNum; index++)
			if (this.views[index] != null)
				sb.append(index).append(':').append(this.views[index].toString()).append('\n');
		return sb.toString();
	}
	
}
