lib = "./bld/libmmf.so"
FT = Float64

function init_soil_depth_clm(nzg::Int)
  slz = fill(0.0, nzg + 1)
  dz = fill(0.0, nzg)
  @ccall lib.__module_rootdepth_MOD_init_soil_depth_clm(
    nzg::Ref{Cint}, slz::Ref{FT}, dz::Ref{FT})::Cvoid
  slz, dz
end

function init_soil_depth(nzg::Int)
  slz = fill(0.0, nzg + 1)
  dz = fill(0.0, nzg)
  @ccall lib.__module_rootdepth_MOD_init_soil_depth(
    nzg::Ref{Cint}, slz::Ref{FT}, dz::Ref{FT})::Cvoid
  slz, dz
end


# 示例调用
slz, dz = init_soil_depth_clm(10)
slz, dz = init_soil_depth(40) # 最多40层
