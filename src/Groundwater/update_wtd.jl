function find_jwt(zwt, slz)
  j = 1
  for k in 2:nzg
    if zwt < slz[k]
      j = k
      break
    end
  end
  return j
end

# qspring = 0.0
function update_wtd(nzg, slz, dz, zwt, qspring, ∑, θ, θ_eq, soiltextures, s_moi_wtd)

  soiltxt = [slz[i] < -0.3 ? soiltextures[1] : soiltextures[2] for i in 1:(nzg+1)]

  # Case 1: Water table going up (totwater > 0)
  if ∑ > 0.0
    if zwt >= slz[1]
      for k in 2:nzg
        zwt < slz[k] && break
      end

      iwtd = k         # z_{i+1}, SM
      kwtd = iwtd - 1  # z_i    , GW & SM 混合层
      nsoil = soiltxt[kwtd]

      # Maximum water that fits in the layer
      _layer = dz[kwtd] * (θ_sat(nsoil) - θ[kwtd])   # 最大补给量, dz应该修改，这里写的不好。

      if ∑ <= _layer # 补给量仅够补给一层
        θ[kwtd] += ∑ / dz[kwtd]
        θ[kwtd] = min(θ[kwtd], θ_sat(nsoil))

        if θ[kwtd] > θ_eq[kwtd]
          sy = θ_sat(nsoil) - θ_eq[kwtd]
          zwt = (θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * slz[iwtd] + θ_sat(nsoil) * slz[kwtd]) / sy
          zwt = min(zwt, slz[iwtd])
        end
        ∑ = 0.0
      else
        # Water enough to saturate the layer
        θ[kwtd] = θ_sat(nsoil)
        ∑ -= _layer
        k1 = iwtd
        for k in k1:(nzg+1)
          zwt = slz[k]
          iwtd = k + 1
          k == nzg + 1 && break

          nsoil = soiltxt[k]
          _layer = dz[k] * (θ_sat(nsoil) - θ[k])
          if ∑ <= _layer
            θ[k] += ∑ / dz[k]
            θ[k] = min(θ[k], θ_sat(nsoil))

            # 这里比较合理的地方，每一层补给过程汇中，地下水的水位是动态上升的。
            if θ[k] > θ_eq[k]
              sy = θ_sat(nsoil) - θ_eq[k]
              zwt = (θ[k] * dz[k] - θ_eq[k] * slz[iwtd] + θ_sat(nsoil) * slz[k]) / sy
              zwt = min(zwt, slz[iwtd])
            end
            ∑ = 0.0
            break
          else
            θ[k] = θ_sat(nsoil)
            ∑ -= _layer
          end
        end
      end
    elseif zwt >= slz[1] - dz[1]
      nsoil = soiltxt[1]
      _layer = (θ_sat(nsoil) - s_moi_wtd) * dz[1]

      if ∑ <= _layer
        s_moi_eq_wtd = θ_sat(nsoil) * (slpots(nsoil) / (slpots(nsoil) - dz[1]))^(1.0 / B(nsoil))
        s_moi_eq_wtd = max(s_moi_eq_wtd, θ_cp(nsoil))

        s_moi_wtd += ∑ / dz[1]
        s_moi_wtd = min(s_moi_wtd, θ_sat(nsoil))
        if s_moi_wtd > s_moi_eq_wtd
          zwt = min((s_moi_wtd * dz[1] - s_moi_eq_wtd * slz[1] + θ_sat(nsoil) * (slz[1] - dz[1])) / (θ_sat(nsoil) - s_moi_eq_wtd), slz[1])
        end
        ∑ = 0.0
      else
        s_moi_wtd = θ_sat(nsoil)
        ∑ -= _layer
        for k in 1:(nzg+1)
          zwt = slz[k]
          iwtd = k + 1
          if k == nzg + 1
            break
          end
          nsoil = soiltxt[k]
          _layer = dz[k] * (θ_sat(nsoil) - θ[k])
          if ∑ <= _layer
            θ[k] = min(θ[k] + ∑ / dz[k], θ_sat(nsoil))
            if θ[k] > θ_eq[k]
              zwt = min((θ[k] * dz[k] - θ_eq[k] * slz[iwtd] + θ_sat(nsoil) * slz[k]) / (θ_sat(nsoil) - θ_eq[k]), slz[iwtd])
            end
            ∑ = 0.0
            break
          else
            θ[k] = θ_sat(nsoil)
            ∑ -= _layer
          end
        end
      end
    else
      nsoil = soiltxt[1]
      _layer = (θ_sat(nsoil) - s_moi_wtd) * (slz[1] - dz[1] - zwt)
      if ∑ <= _layer
        zwt += ∑ / (θ_sat(nsoil) - s_moi_wtd)
        ∑ = 0.0
      else
        ∑ -= _layer
        zwt = slz[1] - dz[1]
        _layer = (θ_sat(nsoil) - s_moi_wtd) * dz[1]
        if ∑ <= _layer
          s_moi_eq_wtd = θ_sat(nsoil) * (slpots(nsoil) / (slpots(nsoil) - dz[1]))^(1.0 / B(nsoil))
          s_moi_eq_wtd = max(s_moi_eq_wtd, θ_cp(nsoil))

          s_moi_wtd += ∑ / dz[1]
          s_moi_wtd = min(s_moi_wtd, θ_sat(nsoil))
          zwt = (s_moi_wtd * dz[1] - s_moi_eq_wtd * slz[1] + θ_sat(nsoil) * (slz[1] - dz[1])) / (θ_sat(nsoil) - s_moi_eq_wtd)
          ∑ = 0.0
        else
          s_moi_wtd = θ_sat(nsoil)
          ∑ -= _layer
          for k in 1:(nzg+1)
            zwt = slz[k]
            iwtd = k + 1
            if k == nzg + 1
              break
            end
            nsoil = soiltxt[k]
            _layer = dz[k] * (θ_sat(nsoil) - θ[k])

            if ∑ <= _layer
              θ[k] += ∑ / dz[k]
              θ[k] = min(θ[k], θ_sat(nsoil))
              if θ[k] > θ_eq[k]
                zwt = (θ[k] * dz[k] - θ_eq[k] * slz[iwtd] + θ_sat(nsoil) * slz[k]) / (θ_sat(nsoil) - θ_eq[k])
              end
              ∑ = 0.0
              break
            else
              θ[k] = θ_sat(nsoil)
              ∑ -= _layer
            end
          end
        end
      end
    end

    # Water springing at the surface
    qspring = ∑

    # Case 2: Water table going down (totwater < 0)
  elseif ∑ < 0.0
    if zwt >= slz[1]
      for k in 2:nzg
        if zwt < slz[k]
          break
        end
      end
      iwtd = k

      k1 = iwtd - 1
      for kwtd in k1:-1:1
        nsoil = soiltxt[kwtd]

        # Max water that the layer can yield
        maxwatdw = dz[kwtd] * (θ[kwtd] - θ_eq[kwtd])

        if -∑ <= maxwatdw
          θ[kwtd] += ∑ / dz[kwtd]
          if θ[kwtd] > θ_eq[kwtd]
            zwt = (θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * slz[iwtd] + θ_sat(nsoil) * slz[kwtd]) / (θ_sat(nsoil) - θ_eq[kwtd])
          else
            zwt = slz[kwtd]
            iwtd -= 1
          end
          ∑ = 0.0
          break
        else
          zwt = slz[kwtd]
          iwtd -= 1
          if maxwatdw >= 0.0
            θ[kwtd] = θ_eq[kwtd]
            ∑ += maxwatdw
          end
        end
      end

      if iwtd == 1 && ∑ < 0.0
        nsoil = soiltxt[1]
        s_moi_eq_wtd = θ_sat(nsoil) * (slpots(nsoil) / (slpots(nsoil) - dz[1]))^(1.0 / B(nsoil))
        s_moi_eq_wtd = max(s_moi_eq_wtd, θ_cp(nsoil))

        maxwatdw = dz[1] * (s_moi_wtd - s_moi_eq_wtd)

        if -∑ <= maxwatdw
          s_moi_wtd += ∑ / dz[1]
          zwt = max((s_moi_wtd * dz[1] - s_moi_eq_wtd * slz[1] + θ_sat(nsoil) * (slz[1] - dz[1])) / (θ_sat(nsoil) - s_moi_eq_wtd), slz[1] - dz[1])
        else
          zwt = slz[1] - dz[1]
          s_moi_wtd += ∑ / dz[1]
          dzup = (s_moi_eq_wtd - s_moi_wtd) * dz[1] / (θ_sat(nsoil) - s_moi_eq_wtd)
          zwt -= dzup
          s_moi_wtd = s_moi_eq_wtd
        end
      end
    elseif zwt >= slz[1] - dz[1]
      nsoil = soiltxt[1]
      s_moi_eq_wtd = θ_sat(nsoil) * (slpots(nsoil) / (slpots(nsoil) - dz[1]))^(1.0 / B(nsoil))
      s_moi_eq_wtd = max(s_moi_eq_wtd, θ_cp(nsoil))

      maxwatdw = dz[1] * (s_moi_wtd - s_moi_eq_wtd)

      if -∑ <= maxwatdw
        s_moi_wtd += ∑ / dz[1]
        zwt = max((s_moi_wtd * dz[1] - s_moi_eq_wtd * slz[1] + θ_sat(nsoil) * (slz[1] - dz[1])) / (θ_sat(nsoil) - s_moi_eq_wtd), slz[1] - dz[1])
      else
        zwt = slz[1] - dz[1]
        s_moi_wtd += ∑ / dz[1]
        dzup = (s_moi_eq_wtd - s_moi_wtd) * dz[1] / (θ_sat(nsoil) - s_moi_eq_wtd)
        zwt -= dzup
        s_moi_wtd = s_moi_eq_wtd
      end
    else
      nsoil = soiltxt[1]
      wgpmid = θ_sat(nsoil) * (slpots(nsoil) / (slpots(nsoil) - (slz[1] - zwt)))^(1.0 / B(nsoil))
      wgpmid = max(wgpmid, θ_cp(nsoil))
      s_yield_dw = θ_sat(nsoil) - wgpmid
      wtdold = zwt
      zwt = wtdold + ∑ / s_yield_dw
      s_moi_wtd = (s_moi_wtd * (slz[1] - wtdold) + wgpmid * (wtdold - zwt)) / (slz[1] - zwt)
    end
    qspring = 0.0
  end

end
