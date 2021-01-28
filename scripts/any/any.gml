///@func any(value,var...)
///@desc Returns true if any variable is equal to the value
///@arg value
///@arg var...
function any() {

	for(var i=1; i<argument_count; i++)
	    if(argument[i] == argument[0])
	        return true;

	return false;


}
