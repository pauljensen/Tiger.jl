% COBRA_MODEL_EXTENDED  Test model in COBRA format

cobra.description = 'COBRA test model (extended)';

%             1  2  3  4  5  6  7  8  9 10
cobra.S =  [  1  0  0  0 -1  0  0  0  0  0;   % A
              0  1  0 -1  0  0  0  0  0  0;   % C
              0  0  1  0  0 -1  0  0  0  0;   % F
              0  0  0  0  1  1 -1  0  0  0;   % B
              0  0  0  1  0  0 -1 -1  0  0;   % D
              0  0  0  0  0  0  1  0 -1  0;   % E
              0  0  0  0  0  0  0  1  0  1 ]; % G

cobra.lb = [ -1 -1 -1 -1 -1  0 -1  0  0 -1 ]';
cobra.ub = [  1  1  1  1  1  1  1  1  1  1 ]';
cobra.rev = cobra.lb < 0;

cobra.c  = [  0  0  0  0  0  0  0  0  1  0 ]';

[m,n] = size(cobra.S);

cobra.b = zeros(m,1);

cobra.rxns = array2names('r%i',1:n);
cobra.rxnNames = cobra.rxns;
cobra.mets = {'A','C','F','B','D','E','G'}';
cobra.metNames = cobra.mets;

cobra.genes = {'g4','g5a','g5b', ...
               'g6a', 'g6b','g7a','g7b','g8a','g8b'}';

cobra.grRules = {'';
                 '';
                 '';
                 'g4';
                 'g5a & g5b';
                 'g6a | g6b';
                 'g7a & g7b';
                 'g8a | g8b';
                 '';
                 ''};

cobra.rules = {'';
               '';
               '';
               'x(1)';
               'x(2) & x(3)';
               'x(4) | x(5)';
               'x(6) & x(7)';
               'x(8) | x(9)';
               '';
               ''};

             
%cobra.rules = convert_grRules(cobra);
cobra.rxnGeneMat = make_rxnGeneMat(cobra);

clear m n
