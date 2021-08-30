using HTTP,JSON,Base64,Sockets
import Plots.savefig

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




token_params = Dict(
    "grant_type" => "password",
    "client_id" => client_id,
    "client_secret" => client_secret,
    "email" => email,
    "password" => password
)

"""
    get_teams()
This function retrieves all teams and their IDs the user is member of.  
"""
function get_all_teams()
    if api_running() == true 
    
        r = HTTP.request("GET",SCINOTE_URL*"/api/v1/teams",header(token_tim))

        body = String(r.body)
        _body = JSON.parse(body)
        n_teams = length(_body["data"])
        teams = Dict(_body["data"][i]["attributes"]["name"] => _body["data"][i]["id"] for i in 1:n_teams)
       
    end
end




"""
    get_projects(team::Int64)
This function retrieves all projects and their IDs from the AK Hasselbrink team.
"""
function get_projects(token_tim,team)
    if api_running() == true 
        resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects",header(token_tim))
        body = String(resp.body)
        data=JSON.parse(body)["data"]
        n_projects = length(data)
        projects  = Dict(data[i]["attributes"]["name"]=> data[i]["id"] for i in  1:n_projects)
    end 
end

"""
    get_experiments(project::Int64,team::Int64)
This function retrieves all experiments and their IDs from the specified project
"""
function get_experiments(token_tim,team,project)
    if api_running() == true 
        resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments",header(token_tim))
        body = String(resp.body)
        data=JSON.parse(body)["data"] 
        n_experiments= length(data)
        experiments = Dict(data[i]["attributes"]["name"]=> data[i]["id"] for i in  1:n_experiments)
    end
end

"""
    get_tasks(team::Int64,project::Int64,experiment::Int64)
This function retrieves all tasks and theis IDs from a specific experiment. 
"""
function get_tasks(token_tim,team,project,experiment;show_archived= false)
    if api_running() == true 
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
end


"""
    get_protocols(team::Int64,project::Int64,experiment::Int64,task::Int64)
This function retrieves all protocols and their IDs from a specific experiment. 
"""
function get_protocols(token_tim,team,project,experiment,task)
    if api_running() == true 
        resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols",header(token_tim))
        body = String(resp.body)
        data=JSON.parse(body)["data"]
        
        n_protocols = length(data)

        protocols = [data[i]["id"] for i in 1:n_protocols]
    end
end


"""
    get_steps(team::Int64,project::Int64,experiment::Int64,task::Int64,protocol::Int64)
This function retrieves all steps and their IDs from a specific protocol. 
"""
function get_steps(token_tim,team,project,experiment,task,protocol)
    if api_running() == true 
        resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols/$(protocol)/steps?page[size]=999",header(token_tim))
        body = String(resp.body)
        data=JSON.parse(body)["data"]

        n_steps = length(data)

        steps = Dict(data[i]["attributes"]["name"] => data[i]["id"] for i in 1:n_steps)
    end
end

#"""
#    get_step_table(team::Int64,project::Int64,experiment::Int64,task::Int64,protocol::Int64,step::Int64)
#This function retrieves the table from specific step. 
#Empty cells will be ignored. Be sure to have a proper 
#"""
#function get_step_table(team::Int64,project::Int64,experiment::Int64,task::Int64,protocol::Int64,step::Int64)
#    if api_running() == true
#        resp=HTTP.request("GET",SCINOTE_URL*"/api/v1/teams/$(team)/projects/$(project)/experiments/$(experiment)/tasks/$(task)/protocols/$(protocol)/steps/$(step)/tables",header(token_tim))
#        body = String(resp.body)
#        data=JSON.parse(JSON.parse(body)["data"][1]["attributes"]["contents"])["data"]
#
#        param_names = [data[i][1] for i in 1:length(data) if data[i][1] !== nothing]
#        param_values = [data[i][2] for i in 1:length(data) if data[i][1] !== nothing]
#        param_units = [data[i][3] for i in 1:length(data) if data[i][1] !== nothing]
#
#            return param_names,param_values,param_units
#    end
#end


"""
    post_plot(plot::Plots.Plot{Plots},experiment::AbstractString; name = "This is a plot")
Post the plot on SciNote as an Attachement in the specified experiment. If you dont specify name the name will be the same as experiment \n
If succesfull the function returns the link to the posted attachment \n
Only works for Experiments in:
    AK Hasselbrink -> FemtoLab -> Spectroelectrochemistry
"""
function post_plot(plot,experiment::AbstractString; name = experiment::AbstractString)


    regex = r"(?:-)((\d|\w)*)"
    task_match = match(regex,experiment)
    task_name = task_match.captures[1]
    if api_running() == true 
        token_tim = token()["access_token"]
        task_id = get_tasks(token_tim,team_id,project_id_femto_lab,experiment_id_SEC)[task_name]
        protocol_id = get_protocols(token_tim,team_id,project_id_femto_lab,experiment_id_SEC,task_id)[1]
        step_id = get_steps(token_tim,team_id,project_id_femto_lab,experiment_id_SEC,task_id,protocol_id)[experiment]

        savefig(plot,"./tmp.png")
        file_data = open("./tmp.png") do io
            base64encode(io)
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
        resp = HTTP.request("POST", SCINOTE_URL*"/api/v1/teams/$team_id/projects/$project_id_femto_lab/experiments/$experiment_id_SEC/tasks/$task_id/protocols/$protocol_id/steps/$step_id/attachments",
            header(token_tim),
            JSON.json(attachment)
        )
        body = String(resp.body)
        data=JSON.parse(body)["data"]
        url = data["attributes"]["file_url"]
   end
end
