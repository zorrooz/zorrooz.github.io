/**
 * Vite dev 插件：数据源变更时重生成内容。
 * 监听两层源：content-src（zh 手写）与 cache/en（英文机器层）。翻译输出只写 cache/en，
 * 不反向写回 content-src，故无「触发→翻译→触发」循环风险；
 * 新增/删除文件需重启 dev server（import.meta.glob 在模块 transform 期解析，否则漏检）。
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

      /** 忽略 Obsidian 配置/回收站目录：workspace.json 频繁写入不应触发重生成风暴 */
      const watcher = watch([contentSrcDir, enSrcDir], {
        ignoreInitial: true,
        ignored: /(^|[/\\])\.(obsidian|trash)([/\\]|$)/,
      })
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
