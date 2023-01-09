using Makie

Makie.@Block FavSelect begin
    @forwarded_layout
    buttons::Vector{Button}
    menu::Menu
    textbox::Textbox
    selection::Observable{Vector{Symbol}}
    @attributes begin
        "Controls if multiple syms can be selected at the same time."
        allowmulti = true
        "The horizontal alignment of the block in its suggested bounding box."
        halign = :center
        "The vertical alignment of the block in its suggested bounding box."
        valign = :center
        "The width setting of the block."
        width = Auto()
        "The height setting of the block."
        height = Auto()
        "Controls if the parent layout can adjust to this block's width"
        tellwidth::Bool = false
        "Controls if the parent layout can adjust to this block's height"
        tellheight::Bool = true
        "The align mode of the block in its parent GridLayout."
        alignmode = Inside()
    end
end

function Makie.initialize_block!(fs::FavSelect, favs, syms)
    syms = syms isa Observable ? syms : Observable(syms)
    setfield!(fs, :selection, Observable([favs[1]]))

    fs.buttons = Button[]
    for (i, fav) in enumerate(favs)
        btn = Button(fs.layout[1,i], label=string(fav))
        push!(fs.buttons, btn)
    end
    fs.menu = Menu(fs.layout[1,length(favs)+1], options=syms)
    fs.textbox = Textbox(fs.layout[1,length(favs)+2], width=150)

    for i in 1:length(favs)
        on(fs.buttons[i].clicks) do n
            _update_selection(fs, favs[i])
        end
    end

    on(fs.menu.selection) do sel
        isnothing(sel) || _update_selection(fs, sel)
    end

    fs.textbox.stored_string.ignore_equal_values = true
    on(fs.textbox.stored_string) do str
        str = replace(str, r"[\s+:]"=>s"") # remove : and whitespace
        if occursin(r"^\s*$", str) # empty
            fs.selection[] = Symbol[]
        elseif !fs.allowmulti[]
            fs.selection[] = [Symbol(str)]
        else
            parts = split(str, ','; keepempty=false)
            fs.selection[] = Symbol.(parts)
        end
    end

    on(fs.selection) do sel
        # handle manu
        sel_from_menu = intersect(sel, syms[])
        if isempty(sel_from_menu)
            fs.menu.prompt[] = "nothing selected"
        else
            fs.menu.prompt[] = repr_selection(sel_from_menu) * " selected"
        end
        fs.menu.i_selected[] = 0

        # buttons
        for (i, fav) in enumerate(favs)
            if fav ∈ sel
                fs.buttons[i].buttoncolor = Makie.COLOR_ACCENT[]
            else
                fs.buttons[i].buttoncolor = RGBf(0.94,0.94,0.94)
            end
        end

        rep = repr_selection(sel)
        fs.textbox.displayed_string[] = isempty(rep) ? " " : rep
        # fs.textbox.stored_string.val = rep
    end

    notify(fs.selection)

    return nothing
end

function _update_selection(fs, sym)
    if fs.allowmulti[] && shift_pressed(fs)
        if sym ∉ fs.selection[]
            push!(fs.selection[], sym)
            notify(fs.selection)
        end
    else
        fs.selection[] = [sym]
    end
end

function repr_selection(sel)
    str = reduce(*, string.(sel) .* ", ")
    first(str, length(str)-2)
end

shift_pressed(x::Makie.Block) = shift_pressed(x.parent)
shift_pressed(x::Figure) = shift_pressed(x.scene)
function shift_pressed(sc::Scene)
    keys = sc.events.keyboardstate
    Makie.Keyboard.left_shift ∈ keys || Makie.Keyboard.right_shift ∈ keys
end
