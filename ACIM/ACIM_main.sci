 // Program to calculate the cumulative effects for assessing the impact of anthopogenic pressures affecting the environmental status
// general workflow:
// 1) getting info about the relationships from the LiACAT tool
// 2) building a structure to pack all info related to the corresponding relationships
// 3) filling the structure with info
// 4) visualising the data with graphs
// 5) building models from the corresponding raw data 
// 6) builing the network based on pathways and cause-impact chains by finding the corresponding sources to the targets
// 7) inserting info about environmental data from monitoring stations or test scenarios
// 8) run the model
// 9) check plausibility and correct errors

// ToDo
// Download the tables comprising info about each relationship of interest from LiACAT
// open the LiACAT-files and use "data in columns" to seperate by tab 
// save the data as csv in a directory called "relationships"
// alternatively, use your own data. Please note that a special structure of the table is required (see example)
// define reference values for all targets
// decide if you want to test a limited number of scenarios or if you want to test a range of scenarios by defining a max and min value for each source.
// Decide what R square value you want to set to filter the models


//////////////////////////////////////////////////////////////////////

clear
clearglobal() 
clc
xdel(winsid())

// ------- (Load) Functions --------
exec("models.sci");
exec("node.sci");


function scenes = GroupByScenes(Results)
listScenes = list();
// group by scenes
for i=1:size(InputsEC,1)// the number of scenes

    CollResultscurrScene = [];
    idxcurrSc = find(Results(:,$)==string(i));
    for idx = idxcurrSc
        CollResultscurrScene =[CollResultscurrScene;Results(idx,:)];
    end
    listScenes($+1)=CollResultscurrScene
end
scenes = listScenes;
endfunction

function [commonTargets1,commonTargets2] = FindPairs(ResPerTarget1,ResPerTarget2)
    disp("find pairs")
    disp("this only works for an input vector with 100 values")

    // find the correct pairs based on the target names
    // this will result in a subset of the data
    targetnames1 = [];
    targetnames2 = [];
    for j = 1:length(ResPerTarget2)
        targetnames2 = [targetnames2;ResPerTarget2(j)(1,3)];
    end
    disp("create a table with the targetnames and reference values based on the following list of targetnames ")
    disp("targetnames")
    commonTargets1 = list();
    commonTargets2 = list();
    for i = 1:length(ResPerTarget1)
        currentTarget = ResPerTarget1(i)(1,3);
        disp(currentTarget)
        idxCurrTarget = find(targetnames2==currentTarget);
        if idxCurrTarget <> [] then
            commonTargets2($+1) = ResPerTarget2(idxCurrTarget);
            commonTargets1($+1) = ResPerTarget1(i); 
            resscenT2 = ResPerTarget2(idxCurrTarget)(100,1);
            resscenT1 = ResPerTarget1(i)(100,1)
            if resscenT1 <>resscenT2
                disp("diff found")
                disp(ResPerTarget1(i)(1,3))
            end
        end
    end
endfunction

function ResultsPerTarget = groupByTargets(Results)
// group by targets
ResultsPerTarget = list();
uniqueTargets = unique(Results(:,2));
for i=1:size(uniqueTargets,1)// the number of targets
    currTarget = uniqueTargets(i);
    CurrResPerTarget = [];
    idxcurrTarg = find(Results(:,2)== currTarget);
    for idx = idxcurrTarg// iterate through all results for this target (through all scenes)
        // add the resulst and the scene identifier
        // collect the result value, target value (:,1),the scene name (:,2),
        // the target name (:,3), and the unit(:,4)
        CurrResPerTarget = [CurrResPerTarget;Results(idx,3),Results(idx,$),Results(idx,2),Results(idx,4)];
    end
    ResultsPerTarget(i) = CurrResPerTarget;
end
endfunction


function GraphsPerTarget2(ResultsPerTargetSinglM,ResultsPerTargetCum)
    disp("grahs per target")
    for i = 1:length(ResultsPerTargetSinglM)
        targetResSinglM = strtod(ResultsPerTargetSinglM(i)(:,1));
        targetResCum = strtod(ResultsPerTargetCum(i)(:,1));
        scenes = strtod(ResultsPerTargetSinglM(i)(:,2));

        figure(i)
        scf(i)
        clf()

        if i == 5 then
            plot2d(scenes,targetResSinglM)
            title(string(ResultsPerTargetSinglM(i)(1,3)),"fontsize",4)
            xlabel("scene number","fontsize",4)
            ylabel("microphytobenthos biomass in ug Chla/ dry wt sediment","fontsize",4)
            //a=gca();//axis
            a=get("current_axes")//get the handle of the newly created axes
            a.axes_visible="on"; // makes the axes visible
            a.font_size=4; //set the tics label font size  
            f=get("current_figure")
            f.children.children(1).children.line_style = 3;
            f.children.children(1).children.thickness = 2;
            f.background=8


            plot2d(scenes,targetResCum)
            title(string(ResultsPerTargetCum(i)(1,3)),"fontsize",4)
            xlabel("scene number","fontsize",4)
            ylabel("microphytobenthos biomass in ug Chla/ dry wt sediment","fontsize",4)
            //a=gca();//axis
            a=get("current_axes")//get the handle of the newly created axes
            a.axes_visible="on"; // makes the axes visible
            a.font_size=5; //set the tics label font size  
            f=get("current_figure")
            f.children.children(1).children.line_style = 8;
            f.children.children(1).children.thickness = 2;
            f.background=8

            xs2png(gcf(),"cumSinglMscenes"+string(ResultsPerTargetCum(i)(1,3))+".png")
        else

            plot2d(scenes,targetResSinglM),
            title(string(ResultsPerTargetSinglM(i)(1,3)),"fontsize",4),
            xlabel("scene number","fontsize",4),
            ylabel(string(ResultsPerTargetSinglM(i)(1,3))+ " in " +string(ResultsPerTargetSinglM(i)(1,4)),"fontsize",4)

            //a=gca();//axis
            a=get("current_axes")//get the handle of the newly created axes
            a.axes_visible="on"; // makes the axes visible
            a.font_size=5; //set the tics label font size  
            f=get("current_figure")
            f.children.children(1).children.line_style = 3;
            f.children.children(1).children.thickness = 2;
            f.background=8


            plot2d(scenes,targetResCum),
            title(string(ResultsPerTargetCum(i)(1,3)),"fontsize",4),
            xlabel("scene number","fontsize",4),
            ylabel(string(ResultsPerTargetCum(i)(1,3))+ " in " + string(ResultsPerTargetCum(i)(1,4)),"fontsize",4)

            a=get("current_axes")//get the handle of the newly created axes
            a.axes_visible="on"; // makes the axes visible
            a.font_size=5; //set the tics label font size  
           // a.tight_limits= "on"; // If this property value is "on" axes adapt to fit exactly with the minima and maxima values of the data_bounds.
            f=get("current_figure")
            f.children.children(1).children.line_style = 8;
            f.children.children(1).children.thickness = 2;
            f.background=8

            xs2png(gcf(),"cumSinglMscenes"+string(ResultsPerTargetCum(i)(1,3))+".png")

    end
    end
    xdel(winsid())  
endfunction

function CorrSinglCumM = correlAndGraph(TargetPairssinglM,TargetPairsCum)
    correlationValues = [];
    targetNamesCorrelations = [];
    for i = 1:length(TargetPairssinglM)
        CorrSinglCumM = correl(strtod(TargetPairsCum(i)(:,1)),strtod(TargetPairssinglM(i)(:,1)));
        
        figure(i)
        plot2d(strtod(TargetPairsCum(i)(:,1)),strtod(TargetPairssinglM(i)(:,1))),,
        title(string(TargetPairsCum(i)(1,3)),"fontsize",4),
        xlabel("including cumulative models","fontsize",4),
        ylabel("excluding cumulative models","fontsize",4)
        a=gca(); 
        f=get("current_figure")
        f.background=8
        f.children.children(1).children.line_style = 8;
        f.children.children(1).children.thickness = 2;
        f.background=8 

        if i==5 then
            plot2d(strtod(TargetPairsCum(i)(:,1)),strtod(TargetPairssinglM(i)(:,1))),,            title("microphytobenthos biomass","fontsize",4),
            xlabel("including cumulative models","fontsize",4),
            ylabel("excluding cumulative models","fontsize",4)
            a=gca(); 
            f=get("current_figure")
            f.background=8
            f.children.children(1).children.line_style = 8;
            f.children.children(1).children.thickness = 2;
            f.background=8 
        end

        xs2png(gcf(),"cumSinglMComp"+string(TargetPairsCum(i)(1,3))+".png")
    correlationValues = [CorrSinglCumM;correlationValues];
    targetNamesCorrelations = [TargetPairssinglM(i)(1,3);targetNamesCorrelations];
    end
    pause
    
    VecCorrelations = [targetNamesCorrelations, string(correlationValues)];
    filenameCorrelCumulaSinglM = "CorrelationSinglMCumula"++".csv"
    // save results as csv file
    Correlations = fullfile(".\CorrelationPearson",filenameCorrelCumulaSinglM)
    csvWrite(VecCorrelations,Correlations,";")
    
endfunction

function OutputTarget = calcResPerTargetGroup(TargetGroupsRefDef, TargetName, ResultValue)
    disp("calcResPerTargetGroup")

    if length(ResultValue)>1 then
        disp("2 result values")
        disp(ResultValue)
    end

    global EnvironNodes
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
    TargetGroups = list();
    // take only every fourth entry to get the group definitions and the
    // members (in the other ones the reference values and the 
    // literature information is stored)
    AllGroupsAndMembers = stripblanks(TargetGroupsRefDef(:,1:4:$));
    // only the first line is needed
    TargetGroupsDef = AllGroupsAndMembers(1,:);
    lengthTargetGroupsDef = size(TargetGroupsDef,2);
    // get the members of each group with the same method
    AllGroupMembers = stripblanks(TargetGroupsRefDef(:,1:4:$));
    // delete the first line of the matrix as here are only the headers
    AllGroupMembers(1,:)=[];
    // get the reference values with the same method but all values in a
    // column
    TargetRefValues = stripblanks(TargetGroupsRefDef(:,2:4:$));
    // delete the first line of the matrix as here are only the headers
    TargetRefValues(1,:)=[];
    // as decimal numbers instead of strings
    TargetRefValues=strtod(TargetRefValues); 

    // organize the info in a list of target groups
    for i = 1:lengthTargetGroupsDef
        DefTargetGroup = TargetGroupsDef(i);
        TargetGroup(i).Def = DefTargetGroup;
        TargetGroup(i).Members = AllGroupMembers(:,i);
        TargetGroup(i).Ref = TargetRefValues(:,i);
        TargetGroups($+1) = TargetGroup(i);
    end
    // group for the case of several results for the same target
    ImpactTargetGroup = list();
    disp("start loops in calcResPerTargetGroup")
   // disp(size(TargetName,1))
    // disp(length(TargetGroups))
    disp(TargetName)

    // get the result value for the target
    for j = 1:size(TargetName,1)
        // search in all groups for the target name
        for k = 1:length(TargetGroups)
            //        disp("k")
            //        disp(k)
            indexGroupOfCurrentTarget = find(stripblanks(TargetGroup(k).Members)==stripblanks(TargetName(j)));
            // if the target name is found in the the group member structure
            if indexGroupOfCurrentTarget <> [] then
                // check if a result value was already provided as input
                // if not, get it by searching through the EnvironNodes
                if ResultValue == [] then
                    disp("ResultValue == []")
                    break
                    if EnvironNodes(j).ValueUpdated == %T then
                        // first add here a 'find' to find only the targets 
                        // of the group of environs, then the next steps can 
                        // follow
                        Idx =  find(list2vec(list2vec(EnvironNodes).Name)==TargetName);
                        if Idx == [] the
                            disp("target not found in EnvironNodes -> check name")
                        end
                        ResultToCurrentTarget = EnvironNodes(Idx).Value;
                        ImpactTarget(j).ModResults($+1)= ResultToCurrentTarget;
                    else
                        disp("no updated value of EnvironNode")
                        break
                    end
                else
                    ResultToCurrentTarget = ResultValue(j);
                end
                // get the Name of the current (k) )target group and define it
                // for the group of the target of interest (current j)      
                ImpactTargetGroup(j).GroupName = TargetGroups(k).Def;
                //             disp(TargetGroups(k).Def)
                //             disp("RefValues")
                // get the corresponding reference value (belonging to the
                // target name)
                ImpactTargetGroup(j).RefValue = TargetGroup(k).Ref(indexGroupOfCurrentTarget);
                // disp(list2vec(ImpactTargetGroup(1).RefValues))
                // calculate the ecological quality ratio by:
                // parameter value/ reference value
                disp("value calculated in OutputTargetGroups")
                disp(ResultToCurrentTarget)
                disp("Refvalue")
                disp(TargetGroup(k).Ref(indexGroupOfCurrentTarget))
                ImpactTargetGroup(j).EQR = ResultToCurrentTarget/TargetGroup(k).Ref(indexGroupOfCurrentTarget)
                if ImpactTargetGroup(j).EQR == %inf  then 
                    disp("check result and reference value - reference value must not be 0")
                    disp(ResultToCurrentTarget)
                end
                if ImpactTargetGroup(j).EQR == -%inf  then 
                    disp("check result and reference value - reference value must not be 0")
                    disp(ResultToCurrentTarget)
                end
                
                ImpactTargetGroup(j).TargetName = TargetGroup(k).Members(indexGroupOfCurrentTarget)
            end // end if case index found
        end // end loop target group
    end // end loop target name

    OutputTarget = ImpactTargetGroup;
    
    // generate a message, if the target name is not in the list of reference values
    if isempty(OutputTarget) == %T then
        disp("no reference value defined for target or target name not in the list of members of any target group - target name:")
        disp(TargetName)
    end

    

    // assessmentValue = aim/status for each member
    // calculate mean impact of each target group
    // display the results
endfunction
// identifies the best Model by finding the smallest AIC 
// (Akaike Informations criterium) of a all models (considers the parameter number)
function filteredModels = filterModelsOnR2adj(models, threshold)
    // skip all models below threshold
    filteredModels = list()
    for j = 1:length(models)
        if models(j).R2adj >= threshold then
            filteredModels($+1) = models(j);
        end
    end
endfunction


function [bestModelAIC, bestModelMinError] = findBestModel(models)
    disp("find best model")

    if (length(models) == 0) then
        bestModelAIC = [];
        bestModelMinError = [];
    else
        bestModelAIC = models(1);
        for i = 2:length(models)
            if (models(i).AIC < bestModelAIC.AIC) then
                bestModelAIC = models(i);
            end
        end
        bestModelMinError = models(1);
        for i = 2:length(models)
            if (models(i).MeanSquareError < bestModelMinError.MeanSquareError) then
                bestModelMinError = models(i);
            end
        end
    end
endfunction

function [bestSingleModelsAIC, bestSingleModelsMinError] = findBestSingleModels(models)
    //disp("function single best models")
    // in the list of filtered models there are likely severals models with only
    // one source but with different source names and equal source names.
    // They need to be grouped so that all models belonging to one source 
    // name are compared with regard to the AIC and MinError

    // identify all the models, which belong to a certain source name

    disp("identify best Single models")

    bestSingleModelsAIC = list();
    bestSingleModelsMinError = list();

    if (length(models) == 0) then
        return;
    end

    sourceModels = list2vec(models);
    ModelsVecSources = list2vec(sourceModels.Sources);
    ModelsVecSourceNames = list2vec(ModelsVecSources.Name);
    SourceNames = unique(ModelsVecSourceNames)
    
    if size(SourceNames,1)> 1 then
        disp("several sources")

    end
    
//    for j =1:length(models)
//        disp(models(j).Sources.Name)
//    end

    for i = 1:size(SourceNames,1)
        // the intesity models belonging to one source name are in 
        // the list of indices intensityModels(indexSourceNames)
        // of these models one best model needs to be selected 
        disp(SourceNames(i))
        disp("new sourceName")
        currentSourceModels = list();
        modelIndicesForSourceName = find(ModelsVecSourceNames == SourceNames(i));
        for idx = modelIndicesForSourceName
            currentSourceModels($+1)=sourceModels(idx);
        end

        // find best model
        [bestSingleModelAIC,bestSingleModelMinError] = findBestModel(currentSourceModels);
        bestSingleModelsAIC($+1) = bestSingleModelAIC;
        bestSingleModelsMinError($+1) = bestSingleModelMinError;
    end

    disp("finished find best single models")
endfunction

function initAllUntrainedModels(overallCombineCount)
    global allUntrainedModels

    disp("initAllUntrainedModels")
    allUntrainedModels = getAllModels(getBaseModels(), overallCombineCount);
endfunction

//optimize the model to get the best params for each of the models
function trainedModels = getTrainedModels(models, sampleData)
    trainedModels = list();
    for model = models
        trainedModels($+1) = trainModel(model, sampleData);
    end
endfunction

function [statsValuesModelNodes,ResultsPerTarget,ResultsTargetGroups,ResultsEQRsPerGroup] = displayEnvironNodes(indices)
    global EnvironNodes
    global ModelNodes
    disp("DisplayEnvironNodes")

    VectEnvironNodeNames = [];
    VectEnvironNodeValues = [];
    VectEnvironNodeUnits = [];
    VectEnvironModelNames = []; 
    VectModelSources = [];
    VectModelMSEs = [];
    VectModelAICs = [];
    VectModelR2s = [];
    VectModelR2adjs = [];
    VectModelParams = [];
    VectModelPublication = [];
    VectModelInputCount = [];

    for i = indices'
        node = EnvironNodes(i);
        //disp(i)
        for modelNodeIndex = node.ModelNodeIndices'
            disp(node.Name)
            disp("Node " + node.Name + ": " + string(node.Value) + " " + node.Unit)
            VectEnvironNodeNames = [VectEnvironNodeNames,node.Name];
            VectEnvironNodeValues = [VectEnvironNodeValues,node.Value];
            VectEnvironNodeUnits = [VectEnvironNodeUnits,node.Unit];
            model = ModelNodes(modelNodeIndex).Model;
            disp("  model type " + model.Name)
            disp(" model params ")
            disp(model.Params')
            disp(["  source names: " string(list2vec(list2vec(model.Sources).Name)')])
            disp("ModelSources and stats per model")
            // put the source names in one single field (but (; seperated))
            VectModelSources = [VectModelSources,strcat(string(list2vec(list2vec(model.Sources).Name)), " and ")];
            VectEnvironModelNames = [VectEnvironModelNames,model.Name];
            VectModelMSEs= [VectModelMSEs,model.MeanSquareError];
            VectModelAICs = [VectModelAICs,model.AIC];
            VectModelR2s = [VectModelR2s,model.R2];
            VectModelR2adjs = [VectModelR2adjs,model.R2adj];
            VectModelParams = [VectModelParams,strcat(string(list2vec(model.Params))," , ")];
            VectModelPublication = [VectModelPublication,model.Publication];
            VectModelInputCount = [VectModelInputCount,model.InputCount];
        end
    end

    HeaderModelNodes = ["TargetName" "modelResult" "unit" "sources" " modelName" "MSE" "AIC" "R2" "R2adj" "InputCount" "Params" "Publication"]; //HeaderModelNodes;
    
    statsValuesModelNodes = [HeaderModelNodes;VectEnvironNodeNames',string(VectEnvironNodeValues)',VectEnvironNodeUnits',VectModelSources',VectEnvironModelNames',string(VectModelMSEs)',string(VectModelAICs)',string(VectModelR2s)',string(VectModelR2adjs)',string(VectModelInputCount)',string(VectModelParams)', string(VectModelPublication)'];


    // collect the results per target
    disp("collect results per target")

    VecEnvironNodes = list2vec(EnvironNodes);
    TargetGroup = [];
    TargetNames = [];
    TargetValue = [];
    TargetUnit = [];
    TargetRef = [];
    TargetEQR = [];
    TargetEQRAggregated = [];

    for j = 1:length(VecEnvironNodes)
        // when the vector of a reference value is empty, the environ should not
        // be a target. Exclude these cases
        if VecEnvironNodes(j).Ref <> [] then
            TargetGroup = [TargetGroup,VecEnvironNodes(j).Group(1)];
            TargetValue = [TargetValue,VecEnvironNodes(j).Value(1)];
            TargetNames = [TargetNames,VecEnvironNodes(j).Name(1)];
            TargetUnit = [TargetUnit,VecEnvironNodes(j).Unit(1)];
            TargetEQRAggregated = [TargetEQRAggregated,VecEnvironNodes(j).EQRAggregated];
            if VecEnvironNodes(j).Unit == [] then
                disp(j)
                break
            end
            TargetRef = [TargetRef,VecEnvironNodes(j).Ref(1)];
            if size(TargetRef,1) > 1 then
                disp("Ref>1")
                disp(j)
                break
            end
            TargetEQR = [TargetEQR,VecEnvironNodes(j).EQR];
        end
    end

    ResultsPerTargetHeader = ["groupname" "target" "value" "unit" "Ref" "EQR"];
    ResultsPerTarget = [ResultsPerTargetHeader;string(TargetGroup)',string(TargetNames)',string(TargetValue)',string(TargetUnit)',string(TargetRef)',string(TargetEQRAggregated)'];

    // collect the results for calculations per target group
    disp("collect results per target group")

    ResultTargetGroup = list();
    ResultEQRsPerGroup = list();
    // get the unique groupnames of the environs
    groupNames = unique(list2vec(VecEnvironNodes.Group));
    // create a vector of the correct length (use only one group definition per environ)
    // group dfinitions would otherwise be listed per model
    // this step is needed to get the correct indices
    TargetGroupsOfEnvirons = [];
    for i = 1:length(VecEnvironNodes.Group)
        TargetGroupsOfEnvirons = [TargetGroupsOfEnvirons,VecEnvironNodes.Group(i)(1)];
    end

    for i = 1:size(groupNames,1)
        // find the indices for the environs for each group name
        //disp(groupNames(i))
        indexEnvironsGroupName = find(TargetGroupsOfEnvirons==groupNames(i));

        GroupNames = [];
        EQRs = [];

        for j = 1:size(indexEnvironsGroupName',1)
            NameEnvironNodesGroup = EnvironNodes(indexEnvironsGroupName(j)(1));
            EQRs = [EQRs,EnvironNodes(indexEnvironsGroupName(j)).EQRAggregated];
            GroupNames = [GroupNames,EnvironNodes(indexEnvironsGroupName(j)).Group];
        end
        ResultTargetGroup(i).GroupNames = GroupNames(1);
        ResultTargetGroup(i).EQRMeans = mean(EQRs);
        ResultEQRsPerGroup(i).GroupNames = GroupNames(1);
        ResultEQRsPerGroup(i).EQRsAggregated = [string(EQRs)]; 
    end
    ResultTargetGroupVecNames = list2vec(list2vec(ResultTargetGroup).GroupNames);
    ResultTargetGroupVecEQRMeans = list2vec(list2vec(ResultTargetGroup).EQRMeans);
    ResultsTargetGroupsHeader = ["groupNames" "EQRs"];

    ResultsTargetGroups = [ResultsTargetGroupsHeader;string(ResultTargetGroupVecNames),string(ResultTargetGroupVecEQRMeans)];
    
    // find the maximum size of the EQRs of a group to be able to form a symetric matrix
    currentMaxSize = 0;
    for i = 1:length(ResultEQRsPerGroup)
        currentSize = size(ResultEQRsPerGroup(i).EQRsAggregated,2);
        if size(ResultEQRsPerGroup(i).EQRsAggregated,2) > currentMaxSize then
            currentMaxSize = currentSize;
        end
    end
    maxSize = currentMaxSize;
    ResultsEQRsPerGroupVec = [];
    // now iterate through the groups and fill create cells with NANs to form a symetric matrix
    for i = 1:length(ResultEQRsPerGroup)
        if size(ResultEQRsPerGroup(i).EQRsAggregated,2) < maxSize then
            sizeCurrentEQRs = size(ResultEQRsPerGroup(i).EQRsAggregated,2)+1;
            for j= sizeCurrentEQRs:maxSize
                ResultEQRsPerGroup(i).EQRsAggregated = [string(ResultEQRsPerGroup(i).EQRsAggregated),"NaN" ];
            end
        end
        ResultsEQRsPerGroupVec = [ResultsEQRsPerGroupVec;list2vec(ResultEQRsPerGroup)(i).EQRsAggregated]
    end
    
    ResultsEQRsPerGroupNames = list2vec(list2vec(ResultEQRsPerGroup).GroupNames);
    ResultsEQRsPerGroup = [ResultsEQRsPerGroupNames,ResultsEQRsPerGroupVec];
    disp("end display environ nodes")

endfunction

function displayGraphForModel(model,ModelType)
    // prepare data for graphs
    // übergeben: bestModels(i).Publication
    // Names Sources, Names Target, Units
    disp("displayGraphForModel")

    dataPointCount = 50;

        if length(model)>1 then
            disp("length model > 1 check")
        end
    
        if length(model.Sources) == 1 then
            NewSourceVector = linspace(min(model.Sources.Data),max(model.Sources.Data),dataPointCount)';
            Result = evalModel(model, NewSourceVector)
            figure()
            scf()
            Sourcelabel = [model.Sources.Name+" in "+model.Sources.Unit];
            Targetlabel = [model.Target.Name+" in "+model.Target.Unit]
            plot2d(model.Sources.Data,model.Target.Data,-3),title(string(model.Publication)+ " R2adj: " + string(model.R2adj)+ " R2: " + string(model.R2) + string(model.Name+string(model.Params)), "fontsize",3),xlabel(Sourcelabel,"fontsize",4),ylabel(Targetlabel,"fontsize",4)
            a=gca(); 
            f=get("current_figure")
            f.background=8 
            plot2d(NewSourceVector,list2vec(Result))
            // save the current figure as png file with the name of the publication
            xs2png(gcf(),string(model.Publication+" "+ModelType+" "+model.Name+".png"))
        end
endfunction

function env = createEnviron(data, unit, name)
    env.Data = strtod(string(data)); // to handle strange behavior of min/max functions
    env.Unit = unit;
    env.Name = name;
endfunction

function Rel = createModelsForRelationship(path, filename)
    disp("csvRead " + path + filename)
    RelAllInfo = csvRead(path + filename,";",".","string",[],[],[],0);
    Reldata = csvRead(path + filename,";",".","double",[],[],[],1);

    Rel.Title=RelAllInfo(2,2)
    Rel.Species=RelAllInfo(2,15)

    Rel.Target = createEnviron(Reldata(:,5),..
    RelAllInfo(2,6),..
    RelAllInfo(2,14));

    Rel.ExpTime = createEnviron(Reldata(:,7),..
    RelAllInfo(2,8),..
    "Time");

    Rel.Source = list();
    Rel.Source(1) = createEnviron(Reldata(:,3),..
    RelAllInfo(2,4),..
    RelAllInfo(2,13));
    Rel.Source(2) = createEnviron(Reldata(:,18),..
    RelAllInfo(2,19),..
    RelAllInfo(2,17));
    Rel.Source(3) = createEnviron(Reldata(:,21),..
    RelAllInfo(2,22),..
    RelAllInfo(2,20));
    Rel.Source(4) = createEnviron(Reldata(:,24),..
    RelAllInfo(2,25),..
    RelAllInfo(2,23));
    Rel.Source(5) = createEnviron(Reldata(:,27),..
    RelAllInfo(2,28),..
    RelAllInfo(2,26));

    Rel.Models = list();
    Rel.AllRelevantSources = list();

    disp(Rel.Title)

    // identify all relevant sources (environs containing different values)
    allSources = lstcat(Rel.Source, list(Rel.ExpTime));
    for source = allSources
        disp(source.Data)
        if ~isnan(min(source.Data)) && min(source.Data)~=max(source.Data) then
            Rel.AllRelevantSources($+1) = source;
        end
    end
    disp("count of relevant sources: " + string(length(Rel.AllRelevantSources)))
    disp(["allRelevantSources: " string(list2vec(list2vec(Rel.AllRelevantSources).Name)')])

    /////---------Models and parameter optimization---------/////
    initAllUntrainedModels(length(Rel.AllRelevantSources));

    disp("Train models")
    for combineCount = 1:length(Rel.AllRelevantSources)
        disp("Input Count: " + string(combineCount))
        sourceCombinations = getCombinations(Rel.AllRelevantSources, combineCount);
        for sources = sourceCombinations
            [allValues, sizes] = list2vec(sources.Data);
            sourceValues = matrix(allValues,sizes(1), size(sizes,1));
            trainedModels = createTrainedModels(sourceValues, Rel.Target.Data, combineCount);
            for trainedModel = trainedModels
                trainedModel.Publication = basename(filename);
                trainedModel.Target = Rel.Target;
                trainedModel.Sources = sources;
                Rel.Models($+1) = trainedModel;
            end
        end
    end
endfunction

function trainedModels = createTrainedModels(sources, targets, combineCount)
    global allUntrainedModels
    trainedModels = getTrainedModels(allUntrainedModels(combineCount), [sources, targets]);
endfunction

function allRels = createAllModels(path)
    disp("-------------Create all models-----------------------")

    csv_files = listfiles(path + "*.csv")
    disp('file count: ' + string(size(csv_files,1)))
    disp(csv_files)

    allRels = list();
    maxCombineCount = 1;
    for i=1:size(csv_files,1)

        [p,n,e] = fileparts(csv_files(i));
        name = n+e;
        if (relDataExists(path, name)) then
            allRels($+1) = loadRel(path, name);
        else
            allRels($+1) = createModelsForRelationship(path, name);
            saveRel(allRels($), path, name);
        end
            
        for model = allRels($).Models
            if (model.InputCount > maxCombineCount) then
                maxCombineCount = model.InputCount;
            end
            // make sure all models have the fileds Bsup and Binf (delete later 
            // as this issue was fixed in the script "models")
            if ~(isfield(model, "Binf")) then
                model.Bsup = [];
                models.Binf = [];
            end
            //showModelAndSources(model);
        end
    end
    // refresh all combined model functions (if models loaded)
    initAllUntrainedModels(maxCombineCount);
endfunction

function [allBestModelsAIC,allBestModelsMinErr, bestSinglModelsAIC,bestSinglModelsMinErr] =identifyAllBestModels(allRels, rSquareThreshold)
    disp("-------------Identify all best models-----------------------")

    allBestModelsAIC = list();
    allBestModelsMinErr = list();
    bestSinglModelsAIC = list();
    bestSinglModelsMinErr = list();
    bestIntensityModelsAIC = list();
    bestIntensityModelsMinErr = list();
    for rel = allRels
        disp ("identify best models for " + rel.Title)
        


        AllfilteredModels = filterModelsOnR2adj(rel.Models, rSquareThreshold);
        // delete all models with time as the only input source 
        disp("filtered models")
        for i = 1:length(AllfilteredModels)
            if ~(isfield(AllfilteredModels(i), "Binf")) then
                AllfilteredModels(i).Bsup = [];
                AllfilteredModels(i).Binf = [];
            end
        end
        
        filteredModels = DeleteSolelyTimeModels(AllfilteredModels);
        if (isempty(filteredModels)) then
            continue;
        end
        

        // find one best model with AIC and ErrorMin
        [bestModelsAIC.BestModel, bestModelsMinErr.BestModel] = findBestModel(filteredModels);

        // create a list with only those filtered models with one source
        // for exposure time
        singleTimeModelsFiltered = list();
        for f = 1:length(AllfilteredModels)
            // if the best model involves several sources then for all 
            // these sources additional single models are needed in case of overlap
            // and for comparison with the simple addition
            if size(AllfilteredModels(f).Sources,2) == 1 then
                if grep(AllfilteredModels(f).Sources(1).Name,"Time") <> []then
                    //disp("single time model")
                    //disp(filteredModels(f).Sources.Name)
                    singleTimeModelsFiltered($+1) = filteredModels(f);
                end
            end
            for i = 1:length(AllfilteredModels(f).Sources) 
                if grep(AllfilteredModels(f).Sources(i).Name,"time") <> []then
                    disp("time")
                end
            end
        end
        
        [bestModelsAIC.BestTimeModel,bestModelsMinErr.BestTimeModel]=findBestModel(singleTimeModelsFiltered)
        singleIntensityModelsFiltered = list();
        // find the best intensity models
        for f = 1:length(filteredModels)
            // if the best model involves several sources then for all 
            // these sources additional single models are needed in case of overlap
            // and for comparison with the simple addition
            if size(filteredModels(f).Sources,2) == 1 then
                singleIntensityModelsFiltered($+1) = filteredModels(f);
            end
        end
        for j=1:length(singleIntensityModelsFiltered)
            disp(singleIntensityModelsFiltered(j).Publication)
            disp(singleIntensityModelsFiltered(j).Sources.Name)
        end
        

        [bestModelsAIC.BestIntensityModels,bestModelsMinErr.BestIntensityModels]=findBestSingleModels(singleIntensityModelsFiltered)
        bestIntensityModelsAIC($+1) = bestModelsAIC.BestIntensityModels;
        bestIntensityModelsMinErr($+1) = bestModelsMinErr.BestIntensityModels;
        //disp("check length single intensity models")

        disp("start displayGraphs")


       // store all the results for this relationship in the list
        allBestModelsAIC($+1) = bestModelsAIC;
        allBestModelsMinErr($+1) = bestModelsMinErr; 
        //bestSinglModelsAIC($+1) = bestModelsAIC.BestIntensityModel;
//        bestSinglModelsMinErr($+1)=bestModelsMinErr.BestIntensityModel;
    end


// flatten the list structure to the correct form
bestSinglModelsAIC = list();
for i = 1:length(bestIntensityModelsAIC)
    if bestIntensityModelsAIC(i)<>[]then
        for j=1:length(bestIntensityModelsAIC(i))
            bestSinglModelsAIC($+1)=bestIntensityModelsAIC(i)(j);
        end
    end
end

bestSinglModelsMinErr = list();
for i = 1:length(bestIntensityModelsMinErr)
    if bestIntensityModelsMinErr(i)<>[]then
        for j=1:length(bestIntensityModelsMinErr(i))
            bestSinglModelsMinErr($+1)=bestIntensityModelsMinErr(i)(j);
        end
    end
end

endfunction

function bestModelsFlatted = flatBestModels(bestModelsPerRel)
    
    bestModelsFlatted = list();

    for bestModels = bestModelsPerRel
        if bestModels.BestModel <> [] then
            bestModelsFlatted($+1) = bestModels.BestModel;
        end
        if bestModels.BestTimeModel <> [] then
            bestModelsFlatted($+1) = bestModels.BestTimeModel;
        end
        for k = 1:length(bestModels.BestIntensityModels)
            if bestModels.BestIntensityModels(k) <> [] then
                bestModelsFlatted($+1) = bestModels.BestIntensityModels(k);
            end
        end
    end
endfunction

function singleModelsFlatted = flatSingleModels(singleModels)
    
    singleModelsFlatted = list();
    
    for singleModel = singleModels
        for k = 1:length(singleModel)
            if singleModel(k) <> [] then
                singleModelsFlatted($+1) = singleModel(k);
            end
        end
    end
endfunction


function [sources, targets] = getAllEnvirons(models)
    // transform three times to get the correct depth
    disp("getAllEnvirons")
    // make sure all models have the fileds Bsup and Binf (delete later 
    // as this issue was fixed in the script "models")
    for i = 1:length(models)
        if ~(isfield(models(i), "Binf")) then
            models(i).Bsup = [];
            models(i).Binf = [];
        end
    end

    vecBestModels = list2vec(models);

    vecBestModelsSources = list2vec(vecBestModels.Sources);
    vecBestModelsSourcesNames = list2vec(vecBestModelsSources.Name);
    sourceNames = unique(vecBestModelsSourcesNames);
    sources = list();
    for sourceName = sourceNames'
        sources($+1) = vecBestModelsSources(vecBestModelsSourcesNames==sourceName)(1);
    end

    vecBestModelsTargets = list2vec(vecBestModels.Target);
    vecBestModelsTargetsNames = list2vec(vecBestModelsTargets.Name);
    targetNames = unique(vecBestModelsTargetsNames);
    targets = list();
    for targetName = targetNames'
        targets($+1) = vecBestModelsTargets(vecBestModelsTargetsNames==targetName)(1);
    end
endfunction

function relevantSourceNames = getAllRelevantSourceNamesFromAllRels(allRels)
    // transform twice to get the correct depth
    NamesAllRelevantSources = list();
    for rel=allRels
        NamesAllRelevantSources($+1) = list2vec(list2vec(rel.AllRelevantSources).Name);
    end
    relevantSourceNames = list2vec(NamesAllRelevantSources);
    relevantSourceNames = unique(relevantSourceNames);
endfunction

function contains = vectorContainsElement(vector, element)
        // if the element is time, it should be ignored if it´s already part
    // of the network (because time as a source only occurs for multidimenional
    // models)
    if string(element) == "Time" then
        contains = %F // (it´s not really false but time should be ignored here)
    else
        contains = length(find(vector==element)) > 0
    end
endfunction

function showModelAndSources(model)
    disp("show model")
    if model == [] then
        disp("[]")
    else
        showModel(model);
        disp(["  source names: " string(list2vec(list2vec(model.Sources).Name)')])
        disp(["  source units: " string(list2vec(list2vec(model.Sources).Unit)')])
    end
endfunction

function filename = createRelFileName(name)
    filename = name + '_rel.data';
endfunction

function exist=relDataExists(path, name)
    exist = isfile(path + createRelFileName(name));
endfunction

function saveRel(relationship, path, name)
    disp("save " + path + name)
    save(path + createRelFileName(name), 'relationship');
endfunction

function relationship = loadRel(path, name)
    disp("load " + path + createRelFileName(name))
    load(path + createRelFileName(name), 'relationship');
endfunction

function SubsetOfModels = DeleteSolelyTimeModels(models)
    SubsetOfModels = list();
    for i = 1:length(models)
        // keep all models that do not have time as the only source
        if models(i).InputCount <> 1 || models(i).Sources(1).Name <> "Time"
            SubsetOfModels($+1) = models(i);
        end
    end
    
endfunction

function [bestModels,TimeModels,SinglModels] = getAllBestModels(path, rSquareThreshold, methodBestModel)
    allRels = createAllModels(path);
    disp("start: get all best Models")

    [allBestModelsAIC,allBestModelsMinErr, bestSinglModelsAIC, bestSinglModelsMinErr] = identifyAllBestModels(allRels, rSquareThreshold);
//    disp("all best models grouped")

    // make all rel models global for debugging purposes
    global BestRelModels

if methodBestModel == "AIC" then
    flattedbestModels = flatBestModels(allBestModelsAIC);
    // delete all models with time as the only source
    bestModels = DeleteSolelyTimeModels(flattedbestModels);
//    flattedIntensityModels = flatSingleModels(list2vec(allBestModelsAIC).BestIntensityModels);
//    intensityModels = DeleteSolelyTimeModels(flattedIntensityModels);
    TimeModels = flatSingleModels(list2vec(allBestModelsAIC).BestTimeModel);
    SinglModels = bestSinglModelsAIC;
end
    if methodBestModel == "MinError" then
        disp("Method model selection: MinError")
        flattedbestModels = flatBestModels(allBestModelsMinErr);
        bestModels = DeleteSolelyTimeModels(flattedbestModels);
        TimeModels = flatSingleModels(list2vec(allBestModelsMinErr).BestTimeModel);
        SinglModels = bestSinglModelsMinErr;
    end

    // calculate percentage of Names of the Best Models compared to all relevant sources
    [modelsSources,modelsTargets] = getAllEnvirons(bestModels);
    allRelevantSources = getAllRelevantSourceNamesFromAllRels(allRels);
    PercentBestModelSources = length(modelsSources)/size(allRelevantSources,1)*100;
    disp("Percentage of the number of sources of the Best Models compared to all relevant sources = " + string(PercentBestModelSources) + "%")
pause
endfunction

function buildNetwork(bestModels)
    disp("Building network")

    // Nodes must be global
    global EnvironNodes
    global ModelNodes
    EnvironNodes = list();
    ModelNodes = list();

    // Create nodes for each environ included in the best models
    disp("create environ nodes")

    [modelsSources,modelsTargets] = getAllEnvirons(bestModels);
    modelsSourceNames = list2vec(list2vec(modelsSources).Name);
    modelsTargetNames = list2vec(list2vec(modelsTargets).Name);

    modelEnvirons = [list2vec(modelsSources); list2vec(modelsTargets)];
    modelEnvironNames = list2vec(list2vec(modelEnvirons).Name);
    allEnvironNames = union(modelsSourceNames,modelsTargetNames);
    for environName = allEnvironNames
        environ = modelEnvirons(modelEnvironNames==environName)(1);
        EnvironNodes($+1) = createEnvironNode(environ.Name, environ.Unit);
    end

    // Create all relevant model nodes
    // preparation: sort models on input count then on sample count descending
    m = [];
    for i = 1:length(bestModels)
        m = [m; bestModels(i).InputCount size(bestModels(i).Sources(1).Data, 1) i];
    end

    m = gsort(m,'lr'); // sorting in lexical decreasing order
    sortedModels = list();
    disp("sorting in lexical decreasing order")
    for i = 1:length(bestModels)
        sortedModels($+1) = bestModels(m(i,3));
    end

    relevantSourceNames = list();

    disp("create model nodes")
    for targetNodeIndex = 1:length(EnvironNodes)
        environNode = EnvironNodes(targetNodeIndex);
        // skip non-target nodes
        if (~vectorContainsElement(modelsTargetNames, environNode.Name)) then
            continue;
        end

        sourceNamesForTarget = list();
        for model=sortedModels
             
            if (model.Target.Name == environNode.Name) then
                // check if any model source already included. skip this model in that case.
                modelSourceAlreadyConsidered = %F;
                for source=model.Sources
                    if (vectorContainsElement(list2vec(sourceNamesForTarget), source.Name)) then
                        modelSourceAlreadyConsidered = %T;
                        break;
                    end
                end
                if (modelSourceAlreadyConsidered) then
                    continue;
                end
                // make sure no models will be included with time as the only input source
                if model.InputCount == 1 && model.Sources.Name == "Time" then
                    continue;
                end

                // create model node
                modelNode = createModelNode(model);
//                disp("modelNode")

                ModelNodes($+1) = modelNode;

                // connect model node with target node
                modelNodeIndex = length(ModelNodes);
                addModelNodeToTargetNode(modelNodeIndex,targetNodeIndex)

                // note model sources for target 
                for source = model.Sources
                    sourceNamesForTarget($+1) = source.Name;
                end
            end
        end
        relevantSourceNames = lstcat(relevantSourceNames, sourceNamesForTarget);
    end

    relevantEnvironNames = union(list2vec(relevantSourceNames), modelsTargetNames);

    // filter not used environs
    relevantEnvironNodes = list();
    for relevantEnvironName = relevantEnvironNames
        environNodeIndex = find(list2vec(list2vec(EnvironNodes).Name)==relevantEnvironName)(1);
        relevantEnvironNodes($+1) = EnvironNodes(environNodeIndex);
    end
    
    EnvironNodes = relevantEnvironNodes;
    // connect source nodes with model nodes
    disp("connect source nodes with model nodes")
    for modelNodeIndex = 1:length(ModelNodes)
        model = ModelNodes(modelNodeIndex).Model
        for source = model.Sources
            sourceNodeIndex = find(list2vec(list2vec(EnvironNodes).Name)==source.Name)(1);
            addSourceNodeToModelNode(sourceNodeIndex,modelNodeIndex);

            allRelevantEnvironNodeIndices($+1) = sourceNodeIndex;
        end
    end
endfunction

function [SolelySourcesIndices, InnerNodesIndices, SolelyTargetsIndices] = getNodesIndexSets()
    //disp("extract the node sets solely sources, inner nodes and solely targets")

    SolelySourcesIndices = [];

    global EnvironNodes

    TargetsIndices = [];
    for i=1:length(EnvironNodes)
        if EnvironNodes(i).ModelNodeIndices == [] then
            SolelySourcesIndices($+1) = i;
        else
            TargetsIndices($+1) = i;
        end
    end

    global ModelNodes
    // sources
    SourcesIndices = [];
    for modelNode=ModelNodes
        SourcesIndices = [SourcesIndices;modelNode.SourceNodeIndices];
    end
    SourcesIndices = unique(SourcesIndices);

    SolelyTargetsIndices = (1:length(EnvironNodes))';
    SolelyTargetsIndices(SourcesIndices) = [];
//    disp("get targets")
//    for i = 1:length(SolelyTargetsIndices)
//        disp(EnvironNodes(i).Name)
//    end
//    pause
        
    InnerNodesIndices = (1:length(EnvironNodes))';
    InnerNodesIndices([SolelyTargetsIndices;SolelySourcesIndices]) = [];
endfunction

function inputValues = retrieveInputParametersFromInputBox()
    global EnvironNodes

    [SolelySourcesIndices, InnerNodesIndices, SolelyTargetsIndices] = getNodesIndexSets();

    SolelySourcesNames=string(list2vec(list2vec(EnvironNodes)(SolelySourcesIndices).Name));
    SolelySourcesValues=string(list2vec(list2vec(EnvironNodes)(SolelySourcesIndices).Value));
    SolelySourcesUnit=string(list2vec(list2vec(EnvironNodes)(SolelySourcesIndices).Unit));

    label = strcat([SolelySourcesNames SolelySourcesUnit], " in ", "c");
    message = 'Please enter here values for each source:';

    InputBoxSources = x_mdialog(message, label, SolelySourcesValues);
    inputValues = strtod(InputBoxSources);  
    
endfunction

function [oInputECValues] = retrieveInputValuesFromFile(InputsEC,InputsECStrings)
 disp("retrieveInputValuesFromFile")

    global EnvironNodes
    //
    [SolelySourcesIndices, InnerNodesIndices, SolelyTargetsIndices] = getNodesIndexSets();

    SolelySourcesNames=string(list2vec(list2vec(EnvironNodes)(SolelySourcesIndices).Name));
    ListoInputECValues = list();
    for j=1:size(InputsEC,1)
        currentInputsECValues = InputsEC(j,:);

        InputsECNames = InputsECStrings(1,:);
        indexSourceNames = [];
        currentOrderedInputECValues = [];
        for i = 1:size(SolelySourcesNames,1)
            //disp(SolelySourcesNames(i))
            indexSourceNamesAdd = find(InputsECNames==SolelySourcesNames(i));
            if indexSourceNamesAdd == [] then 
                disp("source not found in input dataset")
                disp(SolelySourcesNames(i))
                pause
            end

            indexSourceNames = [indexSourceNames,indexSourceNamesAdd];
            addCurrentInputValue = currentInputsECValues(indexSourceNamesAdd);
            currentOrderedInputECValues = [currentOrderedInputECValues,addCurrentInputValue];
            ListoInputECValues(j)= currentOrderedInputECValues';
        end
    end
    
    oInputECValues = [];
    for k=1:length(ListoInputECValues)
        oInputECValues = [oInputECValues,ListoInputECValues(k)]
    end
    
    if length(ListoInputECValues(1)) <> size(SolelySourcesNames,1) then
        disp("missing input values or names not consistent")
        break
    end
endfunction

function injectInputValuesIntoTheNetwork(inputValues, aggregationFunction,AllTargetGroupsRefDef)
    disp("Inject input values and update network.")
   global EnvironNodes

    [SolelySourcesIndices, InnerNodesIndices, SolelyTargetsIndices] = getNodesIndexSets();

    // resetting network (only nodes with values set are valid)
    for i = 1:length(EnvironNodes)
        EnvironNodes(i).ValueUpdated = %F;
    end

    for i = 1:length(SolelySourcesIndices)
        setEnvironNodeValue(SolelySourcesIndices(i), inputValues(i));
    end

    for targetIndex = SolelyTargetsIndices'
        updateEnvironNodeValue(targetIndex, aggregationFunction,AllTargetGroupsRefDef);
    end
endfunction

function result = callModelForNode(nodeIndex)
    global EnvironNodes
    global ModelNodes
    modelNode = ModelNodes(nodeIndex);
    data = [];
    for sourceNodeIndex = modelNode.SourceNodeIndices'
        // ensure source node has an already updated value
        updateEnvironNodeValue(sourceNodeIndex);

        data($+1) = EnvironNodes(sourceNodeIndex).Value;
    end
    result = evalModel(modelNode.Model, data');
endfunction

function updateEnvironNodeValue(nodeIndex, func,AllTargetGroupsRefDef)

    disp("updateEnvironNodeValue")

    global EnvironNodes
    global ModelNodes 
    disp(EnvironNodes(nodeIndex).Name)
    //    disp("EnvironNode")
    //    disp(EnvironNodes(nodeIndex))
    //    disp("ModelNodes")
    //    disp(ModelNodes(nodeIndex))
    if ~EnvironNodes(nodeIndex).ValueUpdated then
        newValues = [];
        //RefValueIndex = find(EnvironNodes(nodeIndex).Name == ) 
        for modelNodeIndex = EnvironNodes(nodeIndex).ModelNodeIndices'
            // new calculated 
            disp("newvalues calc")
            disp(ModelNodes(modelNodeIndex).Model.Target.Name)
            disp(ModelNodes(modelNodeIndex).Model.Sources.Name)

            newvalue = callModelForNode(modelNodeIndex);
            if newvalue == %nan then
                disp(EnvironNodes(nodeIndex).Name)
            //    pause
            end
            newValues($+1) = newvalue;
        end

        DiffEQRValues = [];
        EQRValues = [];
        Log10EQRValues =[];
        RefValues = [];
        GroupDef = [];
        for i = 1:length(newValues)
            RefEQRValueTarget = calcResPerTargetGroup(AllTargetGroupsRefDef, ModelNodes(modelNodeIndex).Model.Target.Name,newValues(i))
            // here it's not a mean but only one EQR value as we had only one input
            disp("output in function updateEnvironNodeValue")
            disp("EQR value of one target based on one model")
            disp(list2vec(RefEQRValueTarget).TargetName)
            disp(list2vec(RefEQRValueTarget).GroupName)
            disp(list2vec(RefEQRValueTarget).EQR)
            // the descrease of the EQR value with respect to the reference value
            DiffEQRvalue = 1-(list2vec(RefEQRValueTarget).EQR);
            DiffEQRValues = [DiffEQRValues,DiffEQRvalue];
            EQRValues = [EQRValues,list2vec(RefEQRValueTarget).EQR];
            Log10EQRValues = [Log10EQRValues,log(list2vec(RefEQRValueTarget).EQR)]
            RefValues = [RefValues,list2vec(RefEQRValueTarget).RefValue];
            GroupDef = [GroupDef,list2vec(RefEQRValueTarget).GroupName];
        end

        aggregatedEQRvalue = 1-func(DiffEQRValues);
        if aggregatedEQRvalue < 0 then
            aggregatedEQRvalue = 0;
        end
        newValue = RefValues(1)*aggregatedEQRvalue;
        disp("single EQR values for one target")
        disp(EQRValues)
        disp("aggregated EQR value")
        disp(aggregatedEQRvalue)
        disp("Ref value")
        disp(RefValues(1))
        disp("new value")
        disp(newValue)  

        // attatch the EQR values, LRRs and group definition to the global EnvironNode
        EnvironNodes(nodeIndex).EQR = EQRValues;
        EnvironNodes(nodeIndex).EQRAggregated = aggregatedEQRvalue;
        EnvironNodes(nodeIndex).Log10EQR = Log10EQRValues;
        EnvironNodes(nodeIndex).Ref = RefValues;
        EnvironNodes(nodeIndex).Group = GroupDef;
        
        setEnvironNodeValue(nodeIndex, newValue)
    end
endfunction

function [statsValuesModelNodes,ResultsPerTarget,ResultsTargetGroups, ResultsEQRsPerGroup] = displayNetwork()
    global EnvironNodes
disp("display network")

    [SolelySourcesIndices, InnerNodesIndices, SolelyTargetsIndices] = getNodesIndexSets();

    disp("Input Nodes:")
    ValuesEnvironInputNodes = displayEnvironNodes(SolelySourcesIndices)

    disp("Inner Nodes:")
    ValuesEnvironInnerNodes = displayEnvironNodes(InnerNodesIndices)

    disp("Output Nodes:") 
    //ValuesEnvironOutputNodes,statsValuesEnvironOutputNodes,ResultsPerTarget
    [statsValuesModelNodes,ResultsPerTarget,ResultsTargetGroups,ResultsEQRsPerGroup] = displayEnvironNodes(SolelyTargetsIndices)

    VecEnvironNodes = list2vec(EnvironNodes);
    SolelyTargetNodes = [];
    for SolelyTargetsIndex = SolelyTargetsIndices
        SolelyTargetNodesAdd = VecEnvironNodes(SolelyTargetsIndex).Name;
        SolelyTargetNodes = [SolelyTargetNodes,SolelyTargetNodesAdd];
    end
    SolelyTargetNodes = list2vec(SolelyTargetNodes);

endfunction

function [statsModelsList,ResultsPerTargetSc,OutputEQRsGroups] = TestScenarios(inputValues, aggregationMethod,AllTargetGroupsRefDef,methodNetworkModels)
    disp("testScenarios")

    if methodNetworkModels == 1 then
        methodNetworkModels = "cum";
    else
        methodNetworkModels = "singlM";
    end

    statsModelsList = list();
    ResultsPerTargetSc = list();
    OutputEQRsGroups = list();

    // prepare vectors to collect the results
    VecResultsPerTarget = [];
    VecStatsModelNodes = [];
    VecResultsPerTargetGr = [];
    VecResultsEQRsPerGroup = [];
    for i = 1:size(inputValues,2)
        injectInputValuesIntoTheNetwork(inputValues(:,i), aggregationMethod,AllTargetGroupsRefDef);
        // max values for the sample points/ German Bight, North Sea
        // sources for values see excel sheet Scenarientest_Berechnungen
        disp("before display network")

        [StatsModelNodes, ResultsPerTarget,ResultsPerTargetGr, ResultsEQRsPerGroup]= displayNetwork()
        disp("network")

        // get the headers
        HeaderResultsPerTarget = [ResultsPerTarget(1,:)];
        HeaderStatsModelNodes = [StatsModelNodes(1,:)];
        HeaderResultsPerTargetGr = [ResultsPerTargetGr(1,:)];

        // exclude the headers from sorting
        StatsModelNodes(1,:) = []; 
        currResultsPerTarget = ResultsPerTarget;
        currResultsPerTarget(1,:)= []; 
        currResultsPerTargetGr = ResultsPerTargetGr;
        currResultsPerTargetGr(1,:)=[];
        currResultsEQRsPerGroup = ResultsEQRsPerGroup;     

        // sort the rows in lexicographic order
        StatsModelNodes = gsort(StatsModelNodes,'lr','i');
        currResultsPerTarget = gsort(currResultsPerTarget,'lr','i');
        disp("add scene idents")

        // add a column the scene identification
        sceneColumn = linspace(i,i,size(currResultsPerTarget,1))';
        sceneIdent = [string(sceneColumn)];
        currResultsPerTarget = [currResultsPerTarget,sceneIdent];

        sceneColumnStats = linspace(i,i,size(StatsModelNodes,1))';
        sceneIdentStats = [string(sceneColumnStats)];
        StatsModelNodes = [StatsModelNodes,sceneIdentStats];

        sceneColumnGr = linspace(i,i,size(currResultsPerTargetGr,1))';
        sceneIdentGr = [string(sceneColumnGr)];
        currResultsPerTargetGr = [currResultsPerTargetGr,sceneIdentGr];

        // collect the results of the different scenes in one vector
        VecResultsPerTarget = [VecResultsPerTarget;currResultsPerTarget];
        VecStatsModelNodes = [VecStatsModelNodes;StatsModelNodes];
        VecResultsPerTargetGr = [VecResultsPerTargetGr;currResultsPerTargetGr];
        disp("before EQRs")
        sortedEQRs = [];
        // sort the EQR values in a differnt way

        for j = 1:size(currResultsEQRsPerGroup(:,1),1) // for the number of EQR groups       
            OneEQRGroup = currResultsEQRsPerGroup(j,:)';
            OneEQRGroup(1)=[];
            for k = 1:size(OneEQRGroup,1)
                sortedEQRs = [sortedEQRs;ResultsEQRsPerGroup(j),OneEQRGroup(k),string(i)];
            end
            //disp("scene"+string(i))
        end

        VecResultsEQRsPerGroup = [VecResultsEQRsPerGroup;sortedEQRs];

        // delete the empty column in the group results
        ColumnToDelete = find(ResultsPerTargetGr(1,:)=="empty");
        ResultsPerTargetGr(:,ColumnToDelete) = [];

        // prepate the output as lists
        OutputEQRsGroups(i)=ResultsPerTargetGr;
        ResultsPerTargetSc(i)= ResultsPerTarget;
       
        disp("before save csv files")        
    end

    // add the headers

    VecResultsPerTarget = [HeaderResultsPerTarget,"scene";VecResultsPerTarget];
    VecStatsModelNodes = [HeaderStatsModelNodes,"scene";VecStatsModelNodes];
     
    HeaderSortedEQRs = ["targetGr" "EQRvalue" "scene"];
    VecResultsEQRsPerGroup =[HeaderSortedEQRs;VecResultsEQRsPerGroup];

    // current filename for results per single target
    filenameResultsPerTarget = "ResultsTargetValues"+"scene"+methodNetworkModels+".csv"
    // save results for target environs
    OuputTargetEnvirons = fullfile(".\ResultsOutputTargets",filenameResultsPerTarget)
    csvWrite(VecResultsPerTarget,OuputTargetEnvirons,";");

    // current filename for statistic values
    filenameResultsTargetStats = "ResultsModelsStats"+methodNetworkModels+".csv"
    OuputTargetsStats = fullfile(".\ResultsOutputTargets",filenameResultsTargetStats)
    csvWrite(VecStatsModelNodes,OuputTargetsStats,";");

    filenameResultsTargetGroups = "ResultsPerGroup"+"scene"+string(i)+methodNetworkModels+".csv"
    OuputGroupsTargets = fullfile(".\ResultsOutputTargets",filenameResultsTargetGroups)
    // delete the column with empty entries and prepare the output

    csvWrite(ResultsPerTargetGr,OuputGroupsTargets,";");  

    // current filename for EQRgroupValues
    filenameEQRgroupAllValues = "ResultsEQRgroupAllValues"+methodNetworkModels+".csv"
    OuputAllEQRValuesgrouped = fullfile(".\ResultsOutputTargets",filenameEQRgroupAllValues)
    csvWrite(VecResultsEQRsPerGroup,OuputAllEQRValuesgrouped,";");

    // Mean MSE
    disp("mean MSE")
    disp("preliminary -- change later to relatative MSE")
    MSEs_AllModels = StatsModelNodes(:,6);
    // delete the header
    MSEs_AllModels(1)= [];
    MSEs_AllModels = strtod(MSEs_AllModels);//convert to decimal
    mean_MSE_AllModels = mean(MSEs_AllModels);
    disp(mean_MSE_AllModels)
    disp("SD of MSEs")
    SD_MSE_AllModels = nanstdev(MSEs_AllModels);
    disp(SD_MSE_AllModels)

    // Mean AIC
    disp("mean AIC")
    disp("preliminary -- change later to relatative AIC")
    AICs_AllModels = StatsModelNodes(:,7);
    AICs_AllModels(1) = [];
    AICs_AllModels = strtod(AICs_AllModels)
    mean_AIC_AllModels = mean(AICs_AllModels);
    disp(mean_MSE_AllModels)
    disp("SD of MSEs")
    SD_AIC_AllModels = nanstdev(AICs_AllModels);
    disp(SD_MSE_AllModels)

    // R2 values
    disp("mean R2")
    R2_AllModels = StatsModelNodes(:,8);
    R2_AllModels(1) = [];
    R2_AllModels = strtod(R2_AllModels);
    mean_R2_AllModels = mean(R2_AllModels);
    disp(mean_R2_AllModels)
    disp("SD of R2")
    SD_R2_AllModels = nanstdev(R2_AllModels);
    disp(SD_R2_AllModels)

    // R2adj values
    disp("mean R2adj")
    R2adj_AllModels = StatsModelNodes(:,9);
    R2adj_AllModels(1) = [];
    R2adj_AllModels = strtod(R2adj_AllModels);
    mean_R2adj_AllModels = mean(R2adj_AllModels);
    disp(mean_R2adj_AllModels)
    disp("SD of R2adj")
    SD_R2adj_AllModels = nanstdev(R2adj_AllModels);
    disp(SD_R2adj_AllModels)

    // put results into lists
    statsModelsList(1).headerMeans = ["meanMSE" "meanAIC" "meanR2" "meanR2adj"];
    statsModelsList(1).means = [mean_MSE_AllModels;mean_AIC_AllModels;...
    mean_R2_AllModels;mean_R2adj_AllModels]';

    statsModelsList(1).headerSDs = ["SDMSE" "SDAIC" "SDR2" "SDR2adj"];
    statsModelsList(1).SDs = [SD_MSE_AllModels;SD_AIC_AllModels;...
    SD_R2_AllModels;SD_R2adj_AllModels]';

    disp("end test scenes")
endfunction

function bestModels = main()
    path = '.\relationships\';
    // get all info from csv files needed
    inputtype = 2; // choose the type of input
    // "1" for three different scenarios and "2" for a range of scenarios
    // based on min and max values 

    if inputtype == 1 then
        // Test scenarios data
        InputsEC = csvRead('.\InputValuesSeagrass\InputValuesTestData.csv',";",".","double",[],[],[],1);
        InputsECStrings = csvRead('.\InputValuesSeagrass\InputValuesTestData.csv',";",".","string",[],[],[],0);
    end
    if inputtype == 2 then
        // Literature min and max values
        InputsECLit = csvRead('.\InputValuesSeagrass\InputValuesTestLitData.csv',";",".","double",[],[],[],1);
        InputsECStrings = csvRead('.\InputValuesSeagrass\InputValuesTestLitData.csv',";",".","string",[],[],[],0);
        InputsECLitnew=[];
        for i=1:size(InputsECLit,2)
            InputsECLitnew= [InputsECLitnew;linspace(InputsECLit(:,i)(1),InputsECLit(:,i)($),100)]; //  points are created
        end
        InputsEC = InputsECLitnew';
    end

    if inputtype == 3 then
        // Reference values ("ideal")
        InputsEC = csvRead('.\InputValuesSeagrass\InputValuesRefData.csv',";",".","double",[],[],[],1);
        InputsECStrings = csvRead('.\InputValuesSeagrass\InputValuesRefData.csv',";",".","string",[],[],[],0);
    end
    TargetGroupsRefDef = csvRead('.\RefValuesAndGroups\AllTargetsGroupsAndRefValues.csv',";",".","string",[],[],[],0);
    disp("start target group")

    // define methods
    rSquareThreshold = 0.6;//currently R2 adjusted value used

    // decide on the method to use to identify the best model
    // the availabel methods are: "AIC" and "MinError" (min least square error)
    methodBestModel = "AIC"; //"MinError";
    aggregationMethod = sum; // sum, prod    

    [allbestModels,bestTimeModels,bestSinglModels] = getAllBestModels(path, rSquareThreshold,methodBestModel);
    
    bestModels = DeleteSolelyTimeModels(allbestModels);

    // add the single models to the best models to fill up gaps of the network
    bestModels = lstcat(bestModels,bestSinglModels);
    
    disp("models ready")


    Models = list(bestModels,bestSinglModels);
    // test model with all cumulative effects (1) and with only single 
    // relationships
    
    ResultsTargets = list();
    ResultsTargetGroups = list();
    for i = 1:2  
        //disp(Models(i))

        buildNetwork(Models(i));
        disp("network ready")

        // Use network
        // define the reference conditions and target groups
        AllTargetGroupsRefDef = csvRead('.\RefValuesAndGroups\AllTargetsGroupsAndRefValues.csv',";",".","string",[],[],[],0);

        // ******* Test scenarios *******
        //inputValues = retrieveInputParametersFromInputBox();

        [oInputECValues] =  retrieveInputValuesFromFile(InputsEC,InputsECStrings)
        inputValues = oInputECValues;

        disp("test scenarios")

        [stats,TargetResults,OutputEQRsGroups] = TestScenarios(inputValues, aggregationMethod,AllTargetGroupsRefDef,i)
        disp("scenarios ready")
        ResultsTargets($+1) = TargetResults;
        ResultsTargetGroups($+1) = OutputEQRsGroups;
        disp("outputs ready for network model")
        if i == 1 then
            disp("end scenarios tests including cumulative models")
        end
        //pause
        if i == 2 then
            disp("end scenarios tests excluding cumulative models")
        end
    end

    xdel(winsid())

    // data analysis- comparison between the methods: cum and singlM
    disp("data analysis difference methods cum singlM")
    //pause
    // **** delete this part later
    InputsECLit = csvRead('.\InputValuesSeagrass\InputValuesTestLitData.csv',";",".","double",[],[],[],1);
    InputsECStrings = csvRead('.\InputValuesSeagrass\InputValuesTestLitData.csv',";",".","string",[],[],[],0);
    InputsECLitnew=[];
    for i=1:size(InputsECLit,2)
        InputsECLitnew= [InputsECLitnew;linspace(InputsECLit(:,i)(1),InputsECLit(:,i)(2),100)]; //  points are created
    end
    InputsEC = InputsECLitnew';
    //*****************************++

    // upload data for cum and singlM and attach a column for the method
    cumResults = csvRead('.\ResultsOutputTargets\ResultsTargetValuesscenecum.csv',";",".","string",[],[],[],1);

    SinglMResults = csvRead('.\ResultsOutputTargets\ResultsTargetValuesscenesinglM.csv',";",".","string",[],[],[],1);

    xdel(winsid())

    // group results
    SinglMScenes = GroupByScenes(SinglMResults);
    singlMResPerTarget = groupByTargets(SinglMResults);
    cumScenes = GroupByScenes(cumResults);
    cumResPerTarget = groupByTargets(cumResults);

    // test if the results correlate with each other (grouped per target)
    // find the correct pairs
    [TargetPairsCum,TargetPairssinglM] = FindPairs(cumResPerTarget,singlMResPerTarget);

    // test the correlation and make corresponding graphs
    CorrSinglCumM = correlAndGraph(TargetPairssinglM,TargetPairsCum)

    GraphsPerTarget2(TargetPairssinglM,TargetPairsCum)

    xdel(winsid())

    disp("summary graph")
    //pause

    // create here a "for loop" for the method comparison
    for k = 1:length(ResultsTargets) // here network build on models from only one input variable (1) and network build including cumulative models (2)
        OutputEQRsGroups = ResultsTargetGroups(k)
        TargetResults = ResultsTargets(k) 

        // for the number of scenes
        for m = 1:length(ResultsTargetGroups(k))
            // figure(k)
            scf(k)
            ResultsPerTargetGr = OutputEQRsGroups(m);
            // delete the empty row
            idx = find(ResultsPerTargetGr(:,1)== "empty");
            ResultsPerTargetGr(idx,:)= [];
            // present the results for changes of the ecosystem elsewhere
            idx = find(ResultsPerTargetGr(:,1)== "biological");
            ResultsPerTargetGr(idx,:)= [];

            // show results in a summary graph
            groupThemes = ResultsPerTargetGr(:,1)
            // delete the first row
            groupThemes(1)=[];
            GroupValues = strtod(ResultsPerTargetGr(:,2))
            GroupValues(1)=[];
            //Bars = strtod(ResultsPerTargetGr(3,:))
            xaxisNumbering = linspace(1,size(GroupValues,1),size(GroupValues,1))'
            //make a scatter plot with filled circles
            if m == 1 then
                scatter(xaxisNumbering,GroupValues,"o")
            end
            if m == 2 then
                scatter(xaxisNumbering,GroupValues,"+")
            end
            if m == 3 then
                scatter(xaxisNumbering,GroupValues,"^")
            end
            ylabel("EQR")
            // add error bars
            //errbar(xaxisNumbering,GroupValues,Bars,Bars)
            if k == 1 then
                xtitle("method 1 - cum")
            else
                xtitle("method 2 - singleM")
            end
            h = gca();
            // keep x-axis flexible but y-axis fixed
            v1 = h.x_ticks.labels;
            GraphXaxisNames = [];

            for i = 1:size(groupThemes,1)
                GraphXaxisNames = [GraphXaxisNames groupThemes(i) " "]
            end
            // delete the last element so that the number of labels is correct
            GraphXaxisNames($) =[];
            h.x_ticks.labels = GraphXaxisNames';
        end
    end
endfunction



// ------- Functions end ----


main();

global BestRelModels
global EnvironNodes
models = flatBestModels(BestRelModels.AIC);
disp("models")
//pause
xdel(winsid())

