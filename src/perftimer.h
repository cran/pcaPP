#pragma once
#include "time.h"

	class CPerfTimer
	{
	public:

		CPerfTimer () : m_start (clock ()) {}

		int GetTimeMS ()
		{
			return FormatTimeMS (clock ()) ;
		}	

		int FormatTimeMS (clock_t end)
		{
			end -= m_start ;

			if (CLOCKS_PER_SEC != 1000)
				end = (clock_t) ((end  * 1000.0 / CLOCKS_PER_SEC)) ;
			return int (end) ;
			
		}
		
	protected:
		clock_t m_start ;
	} ;
