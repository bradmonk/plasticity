/* File parser.mly */
        %token <string> VAL
        %token <string> IDENT
        %token PLUS ARROW 
        %token COMMA
        %token LPAREN RPAREN
        %token EOL
        %token EOF
        %left COMMA
        %left ARROW             /* lowest precedence */
        %left PLUS              /* medium precedence */
        %start reactions        /* the entry point */
        %type <Types.reaction list> reactions
        %%
        reactions:
            reaction_star EOF                   { $1  }
        ;
        
        reaction_star:
          /* empty */                           { [] }
        | reaction reaction_star                { $1::$2 }
        ;
        
        reaction:
          VAL species ARROW species       { ($1, $2, $4)  }
        ;
                         
        species:
        /* empty */                        { [] }
        | each species                     { $1 :: $2 }
        ;
        
        each:
          IDENT LPAREN names RPAREN        { ($1,$3) }
        ;
        
        names:
         /* empty */                      {  [] }
        | VAL                             { [$1] }      
        | IDENT                           { [$1] }
        | VAL COMMA names                 { $1 :: $3 } 
        | IDENT COMMA names               { $1 :: $3 }
        
 /*  I removed this from the definition of the species    */
 /*  | each                             { [$1] } */
     