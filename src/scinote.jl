using HTTP,JSON,Base64,Sockets,Plots,DataFrames
import StringDistances: Levenshtein

const SCINOTE_URL = "http://titus.phchem.uni-due.de:3000"
const client_id = "7vJ_9ypiSlFDLyWYoF_TGpceceDPxXFDnF01Wp3HPMo"
const client_secret = "RJkDV5Lc5ShILMXN2UiTFPnDw-nV1yT3xaPBGXB7w18"
const team_name = "AK Hasselbrink"
const team_id = "2"
const project_id_femto_lab = "9"
const experiment_id_SEC = "53"
const password = "1a2s3d4f"
const email = "tim.laemmerzahl@uni-due.de"
const DNS_uni_due = ip"132.252.3.10"
token_header = ["Content-Type" => "application/json"]
header(token) = Dict(
    "Content-Type" => "application/json",
    "Authorization" => "Bearer $token"
)
token_params = Dict(
    "grant_type" => "password",
    "client_id" => client_id,
    "client_secret" => client_secret,
    "email" => email,
    "password" => password
)

"""
    check_vpn()
Checks whether you are connected to the university network or not.
"""
function check_vpn()
    n = length(Sockets.getalladdrinfo(Sockets.getnameinfo(DNS_uni_due)))
    if n > 1
        true
    else 
        false
    end
end


"""
    api_running()
Returns if the API is running. 
"""
function api_running()
    if check_vpn() == true

        r = HTTP.request("GET",SCINOTE_URL*"/api/health")
        resp = String(r.body)
        if (resp == "RUNNING") == true
            return true
        else
            error("The API seems not to be running. Get Request returns $(resp)")
        end
    else
        error("It seems you are not connected to the Uni-Due Network. Use the VPN. Or Tim fucked up.")
    end
end

"""
token()

Get token with client_id, client_secret, authorization_code, and redirect_uri.
"""
function token()
    if api_running() == true 
        resp = HTTP.request("POST", SCINOTE_URL*"/oauth/token",token_header,
        JSON.json(token_params))
        body = String(resp.body)
        JSON.parse(body)
    end
end

"""
    get_teams(token)
This function retrieves all teams and their IDs the user is member of.  
"""
function get_all_teams(token_tim)
    r = HTTP.request("GET",SCINOTE_URL*"/api/v1/teams",header(token_tim))
    body = String(r.body)
    _body = JSON.parse(body)
    n_teams = length(_body["data"])
    teams = Dict(_body["data"][i]["attributes"]["name"] => _body["data"][i]["id"] for i in 1:n_teams)   
end




"""
    get_projects(token,team::Int64)
This function retrieves all projects and their IDs from the AK Hasselbrink team.
"""
function get_projects(token_tim,team)
    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_projects = length(data)
    projects  = Dict(data[i]["attributes"]["name"]=> data[i]["id"] for i in  1:n_projects)
end

"""
    get_experiments(token,project::Int64,team::Int64)
This function retrieves all experiments and their IDs from the specified project
"""
function get_experiments(token_tim,team,project) 
    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"] 
    n_experiments= length(data)
    experiments = Dict(data[i]["attributes"]["name"]=> data[i]["id"] for i in  1:n_experiments)
end

"""
    get_tasks(token,team::Int64,project::Int64,experiment::Int64)
This function retrieves all tasks and theis IDs from a specific experiment. 
"""
function get_tasks(token_tim,team,project,experiment;show_archived= false)
    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_tasks = length(data)
    if show_archived == false
        tasks = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_tasks if data[i]["attributes"]["archived"] == false)
    else
        tasks = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_tasks)
    end
end


"""
    get_protocols(token,team::Int64,project::Int64,experiment::Int64,task::Int64)
This function retrieves all protocols and their IDs from a specific experiment. 
"""
function get_protocols(token_tim,team,project,experiment,task)
   resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols",header(token_tim))
   body = String(resp.body)
   data=JSON.parse(body)["data"]
   
   n_protocols = length(data)
   protocols = [data[i]["id"] for i in 1:n_protocols]

end


"""
    get_steps(token,team::Int64,project::Int64,experiment::Int64,task::Int64,protocol::Int64)
This function retrieves all steps and their IDs from a specific protocol. 
"""
function get_steps(token_tim,team,project,experiment,task,protocol)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols/$(protocol)/steps?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_steps = length(data)
    steps = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_steps)
end


"""
    get_attachments(token,team::Int64,project::Int64,experiment::Int64,task::Int64,protocol::Int64,step::Int64)
This function retrieves all attachments and their IDs from a specific step. 
"""
function get_attachments(token_tim,team,project,experiment,task,protocol,step)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols/$(protocol)/steps/$(step)/attachments?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_attachments= length(data)
    attachments = Dict(data[i]["attributes"]["file_name"] => data[i]["id"] for i in 1:n_attachments)
end


"""
get_inventories(token_tim,team)

This function retrieves  all active inventories from the specified team.
"""
function get_inventories(token_tim,team)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/inventories?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_inventories = length(data)
    inventories = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_inventories)
end

"""
get_inventoryitems(token_tim,team,inventory)

This function retrieves all inventory items fromthe specified inventory.
"""
function get_inventoryitems(token_tim,team,inventory)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/inventories/$(inventory)/items?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_items = length(data)
    items = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_items)
end

"""
    get_inventorycolumns(token_tim,team,inventory)

This function retrieves all columns from specified inventory.
"""
function get_inventorycolumns(token_tim,team,inventory)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/inventories/$(inventory)/columns?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_columns = length(data)
    columns = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_columns)
end

"""
    get_itemcolumns(token_tim,team,inventory,item)

This function retrieves cells from specified item in inventory.
"""
function get_itemcolumns(token_tim,team,inventory,item)

    resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/inventories/$(inventory)/items/$(item)/cells?page[size]=999",header(token_tim))
    body = String(resp.body)
    data=JSON.parse(body)["data"]
    n_columns = length(data)
    columns = Dict(data[i]["attributes"]["column_id"] => data[i]["attributes"]["value"]["text"] for i in 1:n_columns)
end
"""
    list_inventories()
Lists all inventories in Team AK Hasselbrink.
"""
function list_inventories()
    if api_running() == true
        token_tim = token()["access_token"]
        inventories = get_inventories(token_tim,team_id) |> keys |> collect
        sleep(1)
        return inventories
    end
end

"""
    list_items(;inventory = "FemtoLab samples")
Lists all items in the specified inventory, by default the inventory FemtoLab samples is specified. To see all inventories try the function list_inventories().
With get_sample(sample) you can retrieve all the stored data about the sample from Titus.
"""
function list_items(;inventory = "FemtoLab samples")
    if api_running() == true
        token_tim = token()["access_token"]
        inventories = get_inventories(token_tim,team_id)
        inventory_id = try inventories[inventory] catch end
        if inventory_id === nothing
            k_inventories = inventories |> keys |> collect
            most_prob_idx = [Levenshtein()(inventory,_inventories) for _inventories in k_inventories] |> argmin
            most_prob = k_inventories[most_prob_idx]
            error(""" Inventory "$(inventory)" was not found. Did you mean "$(most_prob)"? Try list_inventories() to see all available inventories on Titus.""")
        end
        items = get_inventoryitems(token_tim,team_id,inventory_id) |> keys |> collect |> sort
        sleep(1)
        return items
    end
end

"""
    get_sample(sample)

Retrieve all data on the sample from Scinote.
"""
function get_sample(sample)
    if api_running() == true
        token_tim = token()["access_token"]
        d_inventories = get_inventories(token_tim,team_id)
        d_items = Dict(k => get_inventoryitems(token_tim,team_id,d_inventories[k]) for k in keys(d_inventories))

        inventory_id = ""
        item_id = ""
        for k in keys(d_inventories)
            if try d_items[k][sample] catch end !== nothing
                inventory_id = d_inventories[k]
                item_id      = d_items[k][sample]
            end
        end

        if isempty(inventory_id) || isempty(item_id)
            error("There are no samples with the name $sample. Please create it.")
        end

        d_header_columns = get_inventorycolumns(token_tim,team_id,inventory_id)
        d_item_columns = get_itemcolumns(token_tim,team_id,inventory_id,item_id)

        ## Construct a DataFrame for the sample. The columns names are fetched from d_header_columns the text is fetched from d_item_columns
        df = DataFrame(:Sample => sample)
        for k in keys(d_header_columns)
           df[!,k] = [try d_item_columns[tryparse(Int,d_header_columns[k])] catch e e="Not specified" end]
        end
        return df  
    end
end




"""
    list_experiments()
Lists all experiment in the group Femto Lab.
"""
function list_experiments()
    if api_running() == true
        token_tim = token()["access_token"]
        experiments = get_experiments(token_tim,team_id,project_id_femto_lab) |> keys |> collect
        sleep(1)
        return experiments
    end
end


"""
    list_tasks(experiment::AbstractString)
Retrieve all task in the given experiment.
e.g. list_steps("Alkanethiols")
"""
function list_tasks(experiment::AbstractString)
    if api_running() == true
        token_tim = token()["access_token"]
        experiments  = get_experiments(token_tim,team_id,project_id_femto_lab)
        experiment_id = try experiments[experiment] catch end
        if experiment_id === nothing
            k_experiments = experiments |> keys |> collect
            most_prob_idx = [Levenshtein()(experiment,_experiment) for _experiment in k_experiments] |> argmin
            most_prob = k_experiments[most_prob_idx]
            error(""" Experiment "$(experiment)" was not found. Did you mean "$(most_prob)"? Try list_experiments() to see all available experiments on Titus.""")
        end
        tasks = get_tasks(token_tim,team_id,project_id_femto_lab,experiment_id) |> keys |> collect |> sort
        sleep(1)
        return tasks
    end
    
end
"""
    list_steps(experiment::AbstractString,task::AbstractString)
Retrieve all steps in the given experiment and task.
e.g. list_steps("Alkanethiols","HDT_01")
"""
function list_steps(experiment::AbstractString,task::AbstractString)
    if api_running() == true
        token_tim = token()["access_token"]
        experiments  = get_experiments(token_tim,team_id,project_id_femto_lab)
        experiment_id = try experiments[experiment] catch end
        if experiment_id === nothing
            k_experiments = experiments |> keys |> collect
            most_prob_idx = [Levenshtein()(experiment,_experiment) for _experiment in k_experiments] |> argmin
            most_prob = k_experiments[most_prob_idx]
            error(""" Experiment "$(experiment)" was not found. Did you mean "$(most_prob)"? Try list_experiments() to see all available experiments on Titus.""")
        end
        tasks = get_tasks(token_tim,team_id,project_id_femto_lab,experiment_id)
        task_id = try tasks[task] catch end
        if task_id === nothing
            k_tasks = tasks |> keys |> collect 
            most_prob_idx = [Levenshtein()(task,_task) for _task in k_tasks] |> argmin
            most_prob = k_tasks[most_prob_idx]
            error(""" Task $(task) was not found. Did you mean "$(most_prob)"? Try list_tasks("$experiment") to see all available tasks in given experiment on Titus.""")
        end
        protocol_id = get_protocols(token_tim,team_id,project_id_femto_lab,experiment_id,task_id)[1]
        steps    = get_steps(token_tim,team_id,project_id_femto_lab,experiment_id,task_id,protocol_id) |> keys |> collect |> sort
        sleep(1)
        return steps
    end

end
"""
    post_step(experiment::AbstractString,task::AbstractString,step::AbstractString,comment::AbstractString)
Post a new step to a task of an experiment. 
e.g. post_step("Alkanethiols","HDT_01","001",get_comment(spectrum))

If this doenst work, make sure you specified the right experiment and task. Use list_experiments() and list_tasks()!!!
"""
function post_step(experiment::AbstractString,task::AbstractString,step::AbstractString,comment::AbstractString)
    if api_running() == true
        token_tim = token()["access_token"]
        experiment_id  = try get_experiments(token_tim,team_id,project_id_femto_lab)[experiment] catch end
        if experiment_id === nothing
            error("$(experiment) was not found. Try list_experiments() to see all available experiments on Titus.")
        end
        task_id = try get_tasks(token_tim,team_id,project_id_femto_lab,experiment_id)[task] catch end
        if task_id === nothing
            error("$(task) was not found. Try list_tasks(experiment) to see all available tasks in given experiment on Titus.")
        end
        protocol_id = get_protocols(token_tim,team_id,project_id_femto_lab,experiment_id,task_id)[1]

        new_step  = Dict(
            "data"=> Dict(
                "attributes"=> Dict(
                    "name" => "$(step)", 
                    "description" => "$(comment)"
                    ),
                "type" => "steps"
            )
        )
        resp = HTTP.request(
            "POST", SCINOTE_URL*"/api/v1/teams/$team_id/projects/$project_id_femto_lab/experiments/$experiment_id/tasks/$task_id/protocols/$protocol_id/steps",
            header(token_tim),
            JSON.json(new_step)
        )
        sleep(1)
    body = String(resp.body)
    data=JSON.parse(body)["data"]["attributes"] |> values |> collect 
    println("The step $(data[1]) was succesfully created.")

    end
end



"""
    post_plot(plot,experiment::AbstractString,task::AbstractString,step::AbstractString; name = experiment::AbstractString )
Post the plot on SciNote as an Attachement in the specified experiment. If you dont specify name the name will be the same as experiment \n
If succesfull the function returns the link to the posted attachment \n

    If this doenst work, make sure you specified the right experiment, task and step. Use list_experiments(),list_tasks() and list_tasks() !!!

e.g. post_plot(plot,"Alkanethiols","HDT_01","001", name = "Awesome Plot" )

"""
function post_plot(plot,experiment::AbstractString,task::AbstractString,step::AbstractString; name = step::AbstractString)

    if api_running() == true
        token_tim = token()["access_token"] 
        experiment_id  = try get_experiments(token_tim,team_id,project_id_femto_lab)[experiment] catch end
        if experiment_id === nothing
            error("$(experiment) was not found. Try list_experiments() to see all available experiments on Titus.")
        end
        task_id = try get_tasks(token_tim,team_id,project_id_femto_lab,experiment_id)[task] catch end
        if task_id === nothing
            error("$(task) was not found. Try list_tasks(experiment) to see all available tasks in given experiment on Titus.")
        end
        protocol_id = get_protocols(token_tim,team_id,project_id_femto_lab,experiment_id,task_id)[1]
        step_id = try get_steps(token_tim,team_id,project_id_femto_lab,experiment_id,task_id,protocol_id)[step] catch end
        if step_id === nothing
            error("$(step) was not found. Try list_steps(experiment,task) to see all available steps in given experiment and task on Titus.")
        end
          
        savefig(plot,"./tmp.png")
        file_data = open("./tmp.png") do io
            Base64.base64encode(io)
        end
        rm("./tmp.png",force=true)
        
        attachment = 
            Dict(
                "data"=> Dict(
                    "attributes"=> Dict(
                        "file_name" => "$name.png", 
                        "file_type" => "image/png",
                         "file_data" => file_data
                        ),
                    "type" => "attachments"
                )
            )
        resp = HTTP.request("POST", SCINOTE_URL*"/api/v1/teams/$team_id/projects/$project_id_femto_lab/experiments/$experiment_id/tasks/$task_id/protocols/$protocol_id/steps/$step_id/attachments",
            header(token_tim),
            JSON.json(attachment)
        )
        sleep(1)
        body = String(resp.body)
        data=JSON.parse(body)["data"]
        url = data["attributes"]["file_url"]
   end
end