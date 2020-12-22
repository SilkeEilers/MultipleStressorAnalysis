
function node=createEnvironNode(name, unit)
    node.Type = "Environ";
    node.Name = name;
    node.Unit = unit;
    node.Value = 0;
    node.ModelNodeIndices = [];
    node.ValueUpdated = %F;
    node.EQR = [];
    node.EQRAggregated = [];
    node.Log10EQR = [];
    node.Group = "empty";
    node.Ref = [];
endfunction

function node=createModelNode(model)
    node.Type = "Model";
    node.Model = model;
    node.SourceNodeIndices = [];
endfunction

function addModelNodeToTargetNode(modelNodeIndex,targetNodeIndex)
    global EnvironNodes
    EnvironNodes(targetNodeIndex).ModelNodeIndices($+1) = modelNodeIndex;
endfunction

function addSourceNodeToModelNode(sourceNodeIndex,modelNodeIndex)
    global ModelNodes
    ModelNodes(modelNodeIndex).SourceNodeIndices($+1) = sourceNodeIndex;
endfunction

function setEnvironNodeValue(nodeIndex, value)
    global EnvironNodes
    EnvironNodes(nodeIndex).Value = value;
    EnvironNodes(nodeIndex).ValueUpdated = %T;
endfunction
