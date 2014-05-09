(* File lexer.mll *)
        {
        open Parser        (* The type token is defined in parser.mli *)
        }
        
        let blank = [' ' '\009' '\012']
        let newline = ('\010' | '\013' | "\013\010")
        let float_literal =
  			('-')? ['0'-'9'] ['0'-'9' '_']*
  			('.' ['0'-'9' '_']* )?
  			(['e' 'E'] ['+' '-']? ['0'-'9'] ['0'-'9' '_']*)?
        let letter = ['A'-'Z' '_' '~' 'a'-'z']
		let digit = ['0'-'9']
		let alphanum = digit | letter
		let ident =  alphanum*   
  			
        rule token = parse
            blank          { token lexbuf }     (* skip blanks *)
          | newline        { token lexbuf }
          | ident as lxm   { IDENT(lxm) }
          | float_literal as lxm { VAL(lxm) }
          | "+"            { PLUS }
          | ","            { COMMA } 
          | "-->"          { ARROW }
          | "("            { LPAREN }
          | ")"            { RPAREN }
          | eof            { EOF }
          