/**
 * Vite dev plugin: regenerate content when the data source layers change.
 * 监听两层源：content-src（zh 手写）与 cache/en（英文机器层）。
 * translation 曾经的循环风险已消除：翻译输出不再回写 content-src，
 * 而是写入 cache/en（重建不会反向写它，故不会形成「触发→翻译→触发」循环）。
 * New/removed files require a server restart (import.meta.glob is resolved
 * at module transform time and would otherwise miss them).
 */

import { watch } from 'chokidar'
import type { Plugin, ViteDevServer } from 'vite'

import { contentSrcDir, enSrcDir } from '../scripts/dataConfig.ts'
import { runAllGenerators } from '../scripts/runAllGenerators.ts'

export function contentDev(): Plugin {
  let timer: ReturnType<typeof setTimeout> | null = null
  return {
    name: 'zorrooz:content-dev',
    apply: 'serve',
    async configureServer(server: ViteDevServer) {
      const run = async (label: string) => {
        try {
          await runAllGenerators()
          server.config.logger.info(`[content] regenerated (${label})`)
          server.ws.send({ type: 'full-reload' })
        } catch (err) {
          server.config.logger.error(`[content] regeneration failed: ${(err as Error).message}`)
        }
      }

      await run('startup')

      const watcher = watch([contentSrcDir, enSrcDir], { ignoreInitial: true })
      watcher.on('all', (event) => {
        if (timer) clearTimeout(timer)
        timer = setTimeout(() => {
          const needsRestart = event === 'add' || event === 'unlink'
          run(needsRestart ? 'restart' : 'change').then(() => {
            if (needsRestart) server.restart()
          })
        }, 300)
      })

      server.httpServer?.on('close', () => watcher.close())
    },
  }
}
