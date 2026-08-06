<template>
  <div class="navigation-tree">
    <div v-for="category in navigationTree" :key="category.name" class="tree-group">
      <div class="tree-label">{{ category.name }}</div>
      <ul class="article-list article-list-root">
        <template v-for="dir in category.children" :key="dir.name">
          <li v-if="dir.files && dir.files.length" class="article-item">
            <div class="directory-node">
              <div class="tree-item tree-item--folder">
                {{ dir.name }}
              </div>
              <ul v-if="dir.files && dir.files.length" class="article-list sub-list files-level">
                <li v-for="file in dir.files" :key="file.path" class="article-item">
                  <router-link
                    :to="toArticle(file.path)"
                    class="tree-item tree-item--child"
                    :class="{ 'tree-item--active': isActive(file.path) }"
                  >
                    {{ file.title }}
                  </router-link>
                </li>
              </ul>
              <ul v-if="dir.children && dir.children.length" class="article-list sub-list">
                <li v-for="sub in dir.children" :key="sub.name" class="article-item">
                  <div class="directory-node">
                    <div class="tree-item tree-item--folder">
                      {{ sub.name }}
                    </div>
                    <ul
                      v-if="sub.files && sub.files.length"
                      class="article-list sub-list files-level"
                    >
                      <li v-for="file in sub.files" :key="file.path" class="article-item">
                        <router-link
                          :to="toArticle(file.path)"
                          class="tree-item tree-item--child"
                          :class="{ 'tree-item--active': isActive(file.path) }"
                        >
                          {{ file.title }}
                        </router-link>
                      </li>
                    </ul>
                  </div>
                </li>
              </ul>
            </div>
          </li>
        </template>
        <li v-for="file in category.files" :key="file.path" class="article-item">
          <router-link
            :to="toArticle(file.path)"
            class="tree-item tree-item--child"
            :class="{ 'tree-item--active': isActive(file.path) }"
          >
            {{ file.title }}
          </router-link>
        </li>
      </ul>
    </div>
  </div>
</template>

<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import { useRoute } from 'vue-router'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { loadCategories } from '@/utils/contentLoader'
import {
  articlePathFromUrl,
  joinRoutePathParam,
  toArticleRoutePath,
  toLocalePath,
} from '@/utils/navigation'
import type { CategoryItem } from '@/types/content'

interface TreeFile {
  title: string
  path: string
}

interface TreeDir {
  name: string
  files: TreeFile[]
  children?: TreeDir[]
}

interface TreeCategory {
  name: string
  files: TreeFile[]
  children: TreeDir[]
}

const route = useRoute()

const navigationTree = ref<TreeCategory[]>([])
const currentPath = ref('')
const { data: categoryData } = useLocalizedContent(() => loadCategories(), [])

function toArticle(path: string) {
  return toLocalePath(toArticleRoutePath(path))
}

function isActive(path: string) {
  const current = currentPath.value.replace(/\.md$/, '')
  const articlePath = path.replace(/\.md$/, '')

  const cleanCurrentPath = current.replace(/-en$/, '')
  const cleanArticlePath = articlePath.replace(/-en$/, '')

  return cleanCurrentPath === cleanArticlePath
}

function buildTree() {
  const path = currentPath.value
  if (!path) {
    navigationTree.value = []
    return
  }
  const segs = path.split('/').filter(Boolean)
  if (segs.length < 2) {
    navigationTree.value = []
    return
  }
  const type = segs[0]
  const group = segs[1]

  let targetItem: CategoryItem | null = null
  outer: for (const section of categoryData.value) {
    for (const item of section.items) {
      if (item.name !== group) continue

      const inArticle = (a: { articleUrl: string }) => a.articleUrl.includes(`/article/${type}/`)
      const hasTypeMatch =
        item.articles?.some(inArticle) ||
        item.categories.some((cat) =>
          cat.articles.some((a) => a.articleUrl.includes(`/article/${type}/${group}/`)),
        )

      if (hasTypeMatch) {
        targetItem = item
        break outer
      }
    }
  }
  if (!targetItem) {
    navigationTree.value = []
    return
  }

  const rootFiles: TreeFile[] = []
  const children: TreeDir[] = []

  const toFile = (title: string, articleUrl: string): TreeFile => ({
    title,
    path: articlePathFromUrl(articleUrl),
  })

  targetItem.articles?.forEach((a) => {
    if (a.articleUrl) rootFiles.push(toFile(a.title, a.articleUrl))
  })

  targetItem.categories.forEach((cat) => {
    const files: TreeFile[] = []
    cat.articles.forEach((a) => {
      if (a.articleUrl) files.push(toFile(a.title, a.articleUrl))
    })
    if (files.length) {
      children.push({
        name: cat.title || cat.key,
        files,
      })
    }
  })

  navigationTree.value = [
    {
      name: targetItem.title || targetItem.name || group,
      files: rootFiles,
      children,
    },
  ]
}

function syncPathAndBuild() {
  currentPath.value = joinRoutePathParam(route.params.path)
  buildTree()
}

watch(() => route.params.path, syncPathAndBuild, { immediate: true })
watch(categoryData, () => buildTree())

onMounted(() => {
  buildTree()
})
</script>

<style scoped>
.navigation-tree {
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

.tree-group {
  margin-bottom: var(--sp-6);
}

.tree-label {
  font-size: 12px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  padding: 0 var(--sp-2);
  margin-bottom: var(--sp-2);
  display: flex;
  align-items: center;
  gap: var(--sp-2);
}

.tree-label::after {
  content: '';
  flex: 1;
  height: 1px;
  background: var(--line);
}

.article-list {
  list-style: none;
  padding-left: 0;
  margin: 0;
}

.article-list-root {
  padding-left: 0;
}

.article-item {
  margin-bottom: 1px;
}

.directory-node {
  padding: 0.1rem 0;
}

.tree-item {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: var(--text-sm);
  color: var(--fg-2);
  padding: 7px 10px;
  border-radius: var(--radius-sm);
  margin-bottom: 1px;
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
  cursor: pointer;
  width: 100%;
  text-align: left;
  border: none;
  background: transparent;
  position: relative;
  text-decoration: none;
  line-height: 1.45;
}

.tree-item:hover {
  background: var(--surface-2);
  color: var(--fg);
}

/* 当前文章 hover 保持蓝色，不落入 hover 的黑色 */
.tree-item--active:hover {
  color: var(--primary);
  background: var(--tint-strong);
}

.tree-item--active {
  color: var(--primary);
  font-weight: 600;
  background: var(--tint);
}

.tree-item--folder {
  font-weight: 600;
  color: var(--fg);
  cursor: default;
}

.tree-item--child {
  padding-left: 10px;
  font-size: var(--text-sm);
}

/* 子目录内的文章（二级及更深）比分类直属文章小一档 */
.files-level .tree-item--child {
  font-size: 13px;
}

.sub-list {
  padding-left: 0;
}
</style>
