% COBRA_MODEL_SMALL  Test model in COBRA format

cobra.description = 'COBRA test model (small)';

%             1  2  3  4  5  6  7  8  9 10
cobra.S =  [  1  1  0 -1;   % A
              0  0  1 -1 ]; % B

cobra.lb = [  0  0  0  0 ]';
cobra.ub = [  1  1  1  1 ]';
cobra.rev = cobra.lb < 0;

cobra.c  = [  0  0  0  1 ]';

[m,n] = size(cobra.S);

cobra.b = zeros(m,1);

cobra.rxns = array2names('r%i',1:n);
cobra.rxnNames = cobra.rxns;
cobra.mets = {'A','B'}';
cobra.metNames = cobra.mets;

cobra.genes = {'g1','g2','g3'}';

cobra.grRules = {'g1';
                 'g2';
                 'g3';
                 ''};

cobra.rules = {'x(1)';
               'x(2)';
               'x(3)';
               ''};
             
%cobra.rules = convert_grRules(cobra);
cobra.rxnGeneMat = make_rxnGeneMat(cobra);

clear m n
