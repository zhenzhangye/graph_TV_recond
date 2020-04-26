/* the break point head file
 *
 * used in message.
 *
 */
#ifndef __BREAKPOINT_H__
#define __BREAKPOINT_H__
class BreakPoint{
public:
	//constructor, deconstructor
	BreakPoint();
	BreakPoint(double lambda, double slope_increment, double shift_increment);

	//operator override
	bool	operator< (const BreakPoint& rhs);
	bool	operator<=(const BreakPoint& rhs);
	bool	operator> (const BreakPoint& rhs);
	bool	operator>=(const BreakPoint& rhs);
	bool	operator==(const BreakPoint& rhs);

	double lambda_;
	double slope_;  	//h(lambda) = slope * lambda + shift;
	double shift_;		//increment for min heap. From left to right.
};

#endif
