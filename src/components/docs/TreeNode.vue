<template>
  <li class="article-item">
    <button
      type="button"
      class="tree-item tree-item--folder"
      :class="[
        depth >= 2 ? `tree-item--l${depth}` : '',
        { 'tree-item--branch-open': isExpanded },
      ]"
      :aria-expanded="isExpanded"
      @click="emit('toggle', dir)"
    >
      <span class="tree-item__text">{{ dir.name }}</span>
      <i
        class="fas fa-chevron-right tree-caret"
        :class="{ 'tree-caret--open': isExpanded }"
        aria-hidden="true"
      ></i>
    </button>

    <ul v-show="isExpanded" class="article-list tree-sublist">
      <li v-for="file in dir.files" :key="file.path" class="article-item">
        <router-link
          :to="toArticle(file.path)"
          class="tree-item"
          :class="[`tree-item--l${depth + 1}`, { 'tree-item--active': isActive(file.path) }]"
        >
          {{ file.title }}
        </router-link>
      </li>

      <TreeNode
        v-for="sub in dir.children"
        :key="sub.id"
        :dir="sub"
        :depth="depth + 1"
        :collapsed-ids="collapsedIds"
        :is-active="isActive"
        :to-article="toArticle"
        @toggle="emit('toggle', $event)"
      />
    </ul>
  </li>
</template>

<script setup lang="ts">
import { computed } from 'vue'

export interface TreeFile {
  title: string
  path: string
}

export interface TreeDir {
  id: string
  name: string
  files: TreeFile[]
  children?: TreeDir[]
}

interface TreeNodeProps {
  dir: TreeDir
  depth: number
  collapsedIds: string[]
  // eslint-disable-next-line no-unused-vars -- 类型签名无函数体，参数名仅作文档
  isActive: (path: string) => boolean
  // eslint-disable-next-line no-unused-vars -- 类型签名无函数体，参数名仅作文档
  toArticle: (path: string) => string
}

const props = defineProps<TreeNodeProps>()

const emit = defineEmits<{ toggle: [dir: TreeDir] }>()

const isExpanded = computed(() => !props.collapsedIds.includes(props.dir.id))
</script>
