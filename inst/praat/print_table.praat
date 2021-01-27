

procedure print_table: .tab, .sep$
	select .tab
	#.tab = selected("Table")
	#.sep$ = ";"
	ncols = Get number of columns
	nrows = Get number of rows
	
	writeInfo: ""
	# Column printing
	for c from 1 to ncols
		cname$ = Get column label: c
		appendInfo: cname$ 
		if c < ncols
			appendInfo: .sep$
		else
			appendInfoLine: ""	
		endif
		
	endfor
	# Data printing
	for r from 1 to nrows
		
		if nrows > 0 
			for c from 1 to ncols
				cname$ = Get column label: c	
				curr$ = Get value: r, cname$
				appendInfo: curr$
				if c < ncols
					appendInfo: .sep$
				else
					appendInfoLine: ""
				endif				
			endfor

		endif

	endfor
endproc

# For testing purposes
#pb = Create formant table (Pols & Van Nierop 1973)
#@print_table: pb , ";"
