/* break point cc file
 *
 * define the structure of a break point.
 *
 */
#include "break_point.h"
BreakPoint::BreakPoint(){
	lambda_	= 0;
	slope_	= 0;
	shift_	= 0;
}

BreakPoint::BreakPoint(double lambda, 
											 double slope, double shift):
											 lambda_(lambda), 
											 slope_(slope), shift_(shift)
											 {}

bool BreakPoint::operator<(const BreakPoint& rhs){
	return lambda_<rhs.lambda_;
}

bool BreakPoint::operator<=(const BreakPoint& rhs){
	return lambda_<=rhs.lambda_;
}

bool BreakPoint::operator>(const BreakPoint& rhs){
	return lambda_>rhs.lambda_;
}

bool BreakPoint::operator>=(const BreakPoint& rhs){
	return lambda_>=rhs.lambda_;
}

bool BreakPoint::operator==(const BreakPoint& rhs){
	return lambda_==rhs.lambda_;
}

