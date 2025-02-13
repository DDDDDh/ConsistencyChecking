package src.cn.edu.nju.moon.consistency.model.observation.constructor;

import src.cn.edu.nju.moon.consistency.model.observation.BasicObservation;

/**
 * @description constructor for {@link BasicObservation}
 * 
 * @author hengxin
 * @date 2012-12-8
 *
 * @see {@link BasicObservation}
 */
public interface IBasicObservationConstructor
{
	public abstract BasicObservation construct();
	
	public abstract String get_ob_id();
}
